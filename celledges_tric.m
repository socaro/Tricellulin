
function [region_jcts,jct,truejcts]=celledges_tric(im_jct,jct,s,edgedetect,h,detectjcts,truejcts,imtric)

tic
Pos=getPosition(h);
im_jct=im_jct(Pos(2):Pos(2)+Pos(4),Pos(1):Pos(1)+Pos(3));
imtric=imtric(Pos(2):Pos(2)+Pos(4),Pos(1):Pos(1)+Pos(3));
if ~isempty(jct)
jct=[jct(:,1)-Pos(1) jct(:,2)-Pos(2)];
end
w=length(im_jct(1,:));
ht=length(im_jct(:,1));
%im_jct = adapthisteq(im_jct);

%im_jct = imadjust(im_jct);
%im_jct1 = medfilt2(im_jct,[round(s/15) round(s/15)]);                        % filter noise

ljct=0;



if edgedetect
    %ed=edge(im_jct,'canny',[],1.7);
    ed=edge(im_jct,'log',[],2);
    ed=bwareaopen(ed,round(s/15));
    ed=imdilate(ed,strel('disk',1));
    ed=bwareaopen(ed,100,4);
    ed=imdilate(ed,strel('disk',2));
    edin=~imfill(~ed,'holes');
    %edin=bwareafilt(~ed,[round((s/5)^2) round((1.7*s)^2)]);
    edin=~bwareaopen(~edin,round((s/3)^2));
    edin=bwmorph(edin,'thin',Inf);
    edin=~bwmorph(edin,'spur',20);
    if detectjcts
        [rj, cj, ~, ~] = findendsjunctions(~edin, 0);
    end
    edin=imerode(edin,strel('disk',1));
else
    edin=imbinarize(im_jct,'adaptive');%,'Sensitivity',0.55);
    edin=bwareaopen(edin,20,4);
    edin=imdilate(edin,strel('disk',2));
    %edin=imfill(~edin,'holes');
%     stats=regionprops(edin,'Eccentricity','PixelIdxList');
%     edinnew=zeros(size(edin));
%     for i=1:length(stats)
%         if stats(i).Eccentricity<0.6
%             edinnew(stats(i).PixelIdxList)=1;
%         end
%     end
    edin=bwareaopen(edin,100);
    edin=imerode(~edin,strel('disk',round(s/20)));
    edin=imfill(edin,'holes');
    edin=bwareaopen(edin,round((s/6)^2),4);
    edin=imclose(edin,strel('disk',round(s/20)));
    %edin=bwareafilt(edin,[round((s/5)^2) round((1.7*s)^2)]);
    edin=bwmorph(~edin,'thin',Inf);
    edin=~bwmorph(edin,'spur',20);
    if detectjcts
        [rj, cj, ~, ~] = findendsjunctions(~edin, 0);
    end
    edin=imerode(edin,strel('disk',1));
    %edin=bwareafilt(edin,[round((s/5)^2) round((2*s)^2)]);
    %edin=imclose(edin,strel('disk',round(s/8)));
end

if detectjcts
delete=false(size(rj));

for i=1:length(rj)-1
   [IDX,d] = knnsearch([rj(1:end-i-1),cj(1:end-i-1)],[rj(end-i),cj(end-i)],'K',5); 
   d_idx=find(d<8);
   delete(IDX(d_idx))=true;
   
end

rj=rj(~delete);
cj=cj(~delete);
end
toc

%[L,Num]=bwlabel(edin);
stats=regionprops(edin,'Area','Perimeter','PixelIdxList','Centroid','Extrema');
region=cell(length(stats));
% [region_row,region_col]=find(L);
% region_coord=[region_row region_col];
% [idx,dist]=knnsearch(region_coord,jct,'K',s^2);
% %junction=L(region_row(idx),region_col(idx));
% junction=zeros(size(idx));
% for i=1:length(idx(:,1))
%    junction(i,:)=L(region_row(idx(i,:)),region_col(idx(i,:)));
% end
toc
totalbw=zeros(size(edin));
%shapefactor=NaN(size(edin));
id=0;

% for i=1:length(stats)
%     %if (stats(i).Centroid(1)<length(edin(1,:))-s)&&(stats(i).Centroid(1)>s)&&(stats(i).Centroid(2)<length(edin(:,1))-s)&&(stats(i).Centroid(2)>s)
%     
%     %in=inpolygon(stats(i).Centroid(1),stats(i).Centroid(2),[Pos(1)-s Pos(1)-s Pos(1)+Pos(3)+s Pos(1)+Pos(3)+s Pos(1)-s],[Pos(2)-s Pos(2)+Pos(4)+s Pos(2)+Pos(4)+s Pos(2)-s Pos(2)-s]);
%     %in=inpolygon(stats(i).Centroid(1),stats(i).Centroid(2),[Pos(1) Pos(1) Pos(1)+Pos(3) Pos(1)+Pos(3) Pos(1)],[Pos(2) Pos(2)+Pos(4) Pos(2)+Pos(4) Pos(2) Pos(2)]);
%     %in=inpolygon(stats(i).Extrema(:,1),stats(i).Extrema(:,2),[Pos(1) Pos(1) Pos(1)+Pos(3) Pos(1)+Pos(3) Pos(1)],[Pos(2) Pos(2)+Pos(4) Pos(2)+Pos(4) Pos(2) Pos(2)]);
%     in=inpolygon(stats(i).Extrema(:,1),stats(i).Extrema(:,2),[s/2 s/2 w-s/2 w-s/2 s/2],[s/2 ht-s/2 ht-s/2 s/2 s/2]);
% 
%     if sum(in)==length(in)
% %         sf=stats(i).Perimeter./sqrt(stats(i).Area);
% 
%         currentreg=zeros(size(edin));
%         currentreg(stats(i).PixelIdxList)=1;
%         currentreg=imclose(currentreg,strel('disk',round(s/4)));
%         totalbw(currentreg==1)=1;
%         id=id+1;
%         reg{id}=currentreg;
%         
%         
%         
% %         currentreg=imerode(currentreg,strel('disk',4));
% %         newstats=regionprops(logical(currentreg),'PixelIdxList');
% %         for j=1:length(newstats)
% %             id=id+1;
% %             reg{id}=zeros(size(edin));
% %             reg{id}(newstats(j).PixelIdxList)=1;
% %             reg{id}=imclose(reg{id},strel('disk',round(s/5)));
% %             reg{id}=imdilate(reg{id},strel('disk',1));
% %             %reg{id}=bwmorph(reg{id},'thicken',4);
% %             totalbw(reg{id}==1)=1;
% %             %for imaging shape factor (usually commented)
% %             currentstats=regionprops(reg{id},'Area','Perimeter');
% %             if ~isempty(currentstats)&&stats(i).Area<5*s^2
% %             sf=currentstats.Perimeter./sqrt(currentstats.Area);
% %             if sf>3&&currentstats.Area>(s/4)^2
% %             shapefactor(reg{id}==1)=sf;
% %             end
% %             end
% %         end
%     else
%         totalbw(stats(i).PixelIdxList)=1;
%     end
% end
toc
if detectjcts
%  bwthin=bwmorph(~totalbw,'thin',Inf);
%  [rj, cj, ~, ~] = findendsjunctions(bwthin, 0);

ss=7;
cn=[];
rn=[];
    for i=1:length(cj)
        if cj(i)>ss && cj(i)<w-ss && rj(i)>ss && rj(i)<ht-ss
        imnodes=imtric(rj(i)-ss:rj(i)+ss,cj(i)-ss:cj(i)+ss);
        [cc, rr]=nodesmax(imnodes,2,2);
        if isempty(cc)
            cn=[cn; cj(i)];
            rn=[rn; rj(i)];
        else
            try
                [~,ind]=max(imnodes(sub2ind(size(imnodes),floor(rr),floor(cc))));
            catch
                [~,ind]=max(imnodes(sub2ind(size(imnodes),ceil(rr),ceil(cc))));
            end
            cn=[cn; cc(ind)+cj(i)-ss-1];
            rn=[rn; rr(ind)+rj(i)-ss-1];
        end
        else
            cn=[cn;cj(i)];
            rn=[rn;rj(i)];
        end
    end
cj=cn;
rj=rn;

%injct=inpolygon(cj,rj,[Pos(1)-s Pos(1)-s Pos(1)+Pos(3)+s Pos(1)+Pos(3)+s Pos(1)-s],[Pos(2)-s Pos(2)+Pos(4)+s Pos(2)+Pos(4)+s Pos(2)-s Pos(2)-s]);
if ~isempty(jct)

jctorig=jct;
ljct=length(jctorig(:,1));
%jct=[cj(injct), rj(injct)];
jct=[cj rj];

%truejcts=[truejcts;true(length(jct(:,1)),1)];
else
    jctorig=[];
    %jct=[cj(injct), rj(injct)];
    jct=[cj rj];
    ljct=0;
    truejcts=[];
end
%stats=regionprops(totalbw,'Area','Perimeter','PixelIdxList','Centroid');
end
toc
id=0;
pjcts=[];
for i=1:length(stats)
      in=inpolygon(stats(i).Extrema(:,1),stats(i).Extrema(:,2),[s/2 s/2 w-s/2 w-s/2 s/2],[s/2 ht-s/2 ht-s/2 s/2 s/2]);
      if sum(in)==length(in)&&stats(i).Area>(s/4)^2&&stats(i).Area<(1.7*s)^2
%     in=inpolygon(stats(i).Centroid(1),stats(i).Centroid(2),[Pos(1)+s/2 Pos(1)+s/2 Pos(1)+Pos(3)-s/2 Pos(1)+Pos(3)-s/2 Pos(1)+s/2],[Pos(2)+s/2 Pos(2)+Pos(4)-s/2 Pos(2)+Pos(4)-s/2 Pos(2)+s/2 Pos(2)+s/2]);
%     if in==1
%         [region_row, region_col]=find(reg{i});
%         region_coord=[region_col region_row];
        [region_row, region_col]=ind2sub(size(im_jct),stats(i).PixelIdxList);
        region_coord=[region_col region_row];
%         currentreg=zeros(size(edin));
%         currentreg(stats(i).PixelIdxList)=1;
            id=id+1;
            [idx,dist]=knnsearch(jct,region_coord,'K',2);
            region_jcts_this=unique(idx(dist<s/8));
            if detectjcts
            region_jcts{id}=[];
            for k=1:length(region_jcts_this)
                jctindex=find(pjcts==region_jcts_this(k));
                if ~isempty(jctindex)
                    region_jcts{id}=[region_jcts{id}; jctindex+ljct];
                else
                    region_jcts{id}=[region_jcts{id}; length(pjcts)+1+ljct];
                    pjcts=[pjcts;region_jcts_this(k)];
                end
            end
            else
                region_jcts{id}=region_jcts_this;
            end
%         end
%          
%     end
      end
end
toc
if detectjcts
   %pjcts=unique(pjcts);
   jct=[jctorig;jct(pjcts,1) jct(pjcts,2)]; 
   truejcts=[truejcts;true(length(pjcts),1)];
end

jct=[jct(:,1)+Pos(1) jct(:,2)+Pos(2)];
% for i=1:Num
%         [region_row, region_col]=find(L==i);
%         region_coord=[region_row region_col];
%         region{i}=zeros(size(L));
%         region{i}(region_row,region_col)=1;
% 
%             for j=1:length(jct(:,1))
%                 [idx,dist]=knnsearch(region_coord,jct(j,:));
%                 junction(j,i)=dist;
%             end
% %             if isempty(find(region_coord==1,1))&&isempty(find(region_coord==size(im_jct,1),1))&&isempty(find(region_coord==size(im_jct,2),1))
% %             fullcell(i)=true;
% %             end
% 
%  
% end

% [B_jct,I_jct]=sort(junction,2);
% I_jct_corr=I_jct(:,1:7);
% I_jct_corr(B_jct(:,1:7)>round(s/2))=0;
% 
% for i=1:length(L)
%     [region_jcts{i},~]=find(I_jct_corr==i);
% end


end