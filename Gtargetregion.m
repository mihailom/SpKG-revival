function [deltGtfv] = Gtargetregion(intervals,structure,l)

%Turner 2004 stacking parameters 
table4 = [0,0,0,-0.9;0,0,-2.2,0;0,-2.1,0,-0.6;-1.1,0,-1.4,0];
table7 = [0,0,0,-2.1;0,0,-3.3,0;0,-2.4,0,-1.4;-2.1,0,-2.1,0];
table10 = [0,0,0,-2.4;0,0,-3.4,0;0,-3.3,0,-1.5;-2.2,0,-2.5,0];
table12 = [0,0,0,-1.3;0,0,-2.5,0;0,-2.1,0,-0.5;-1.4,0,+1.3,0];
table13 = [0,0,0,-1.3;0,0,-2.4,0;0,-2.1,0,-1.0;-0.9,0,-1.3,0];
table15 = [0,0,0,-1.0;0,0,-1.5,0;0,-1.4,0,0.3;-0.6,0,-0.5,0];
au_penalty = 0.45;
gguc_energy = -4.12;
%deltGtfv = zeros(length(intervals),1); 
deltGtfv = zeros(size(intervals,1),1);  
for k=1:size(intervals,1)
    i = intervals(k,1);%-21;
    j = intervals(k,2);%-21;
    deltGtf = 0;
    
    if i == 1
        i = 2;
    end
    if j == l
        j = l-1;
    end
    sys = 0;
    m = 0;
    for n = i:(j+1)
    if m == 2; 
        sys = 0;
    end
    if sys == 1;
       m = m + 1;
      continue; 
    end
    if structure(n-1,2)==0
         deltGtf = deltGtf + 0;
         continue;
    elseif structure(n,2)==0 
        if (structure(n-1,1)*2) + structure(structure(n-1,2),1) == 6 || (structure(n-1,1)*2) + structure(structure(n-1,2),1) == 9 ||(structure(n-1,1)*2) + structure(structure(n-1,2),1) == 10 || (structure(n-1,1)*2) + structure(structure(n-1,2),1) == 19
            deltGtf = deltGtf + au_penalty;
        end
    end
    if structure(n,2)==0
         deltGtf = deltGtf + 0;
         continue;
    end
%     if (structure(n-1,1)==3)&&(structure(n,1)==3)&&(structure(n+1,1)==4)&&(structure(n+2,1)==2)
%            if (structure(structure(n-1,2),1)==2)&&(structure(structure(n,2),1)==2)&&(structure(structure(n+1,2),1)==2)&&(structure(structure(n+2,2),1)==2)
%                deltGtf = deltGtf + gguc_energy;
%                 sys = 1;
%                 continue; 
%             end
%      end
            tables = (structure(n-1,1)*2) + structure(structure(n-1,2),1);
        switch tables
            case 6
                if n-1 == 1 || n-1 == i-1
                    deltGtf = deltGtf + au_penalty;
                elseif structure(n-2,2) == 0
                    deltGtf = deltGtf + au_penalty; 
                end
                deltGtf = deltGtf + table4(structure(n,1),structure(structure(n,2),1));
            case 9   
                 if n-1 == 1 || n-1 == i-1
                    deltGtf = deltGtf + au_penalty;
                elseif structure(n-2,2) == 0
                    deltGtf = deltGtf + au_penalty; 
                end
                deltGtf = deltGtf + table13(structure(n,1),structure(structure(n,2),1));
            case 7
                deltGtf = deltGtf + table7(structure(n,1),structure(structure(n,2),1));
            case 8
                deltGtf = deltGtf + table10(structure(n,1),structure(structure(n,2),1));
            case 10    
                 if n-1 == 1 || n-1 == i-1
                    deltGtf = deltGtf + au_penalty;
                elseif structure(n-2,2) == 0
                    deltGtf = deltGtf + au_penalty; 
                end
                deltGtf = deltGtf + table12(structure(n,1),structure(structure(n,2),1));
            case 19
                if n-1 == 1 || n-1 == i-1
                    deltGtf = deltGtf + au_penalty;
                elseif structure(n-2,2) == 0
                    deltGtf = deltGtf + au_penalty; 
                end
                deltGtf = deltGtf + table15(structure(n,1),structure(structure(n,2),1));      
         end
    end
    deltGtfv(k,1)= deltGtf; 
end
