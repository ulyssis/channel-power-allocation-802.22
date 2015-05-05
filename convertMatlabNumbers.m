function [array] = convertMatlabNumbers(matrix)

array = mat2str(matrix, 8);
% disp(array);

% delete the date between 0 and ';', which is (;, 0]
    i=2;
    while (i~=size(array,2))    
        while(array(i)~=';' && i<size(array,2))
            i=i+1;
        end
        
        j=i;
                
        while ~(array(j)=='0'&&array(j-1)==' '&&(array(j+1)==' '||array(j+1)==']')) && (j~= size(array,2))
            j=j+1;
        end
        
        for k=i:j
            array(k)=' ';
        end

    end


array(1)='';
array(1)='';
array(size(array,2))='';

% add comma between numbers
i=1;
while(array(i)~='.' && i<size(array,2))
    i=i+1;
end
for k=i+1:1:size(array,2)
    if(array(k)=='.')
       m=0;
       while(array(k-m)~=' ')
          m=m+1; 
       end
       array(k-m)=',';
    end
end
