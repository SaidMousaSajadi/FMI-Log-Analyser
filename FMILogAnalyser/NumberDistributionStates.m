function CalU = NumberDistributionStates(Base,Target)
    % Example:
    % Base = 5 ; % length of rows: The first array has the coefficient of 1 and the last array of each line has the Base coefficient.
    % Target = 5 ; % Target is the number you want to write as Base vectors.
    % Output = NumberDistributionStates(Base,Target)
    % >> Output :
    % [5     0     0     0     0 ; ... *1     *2     *3     *4     *5  :  Base
    %  3     1     0     0     0 ;
    %  2     0     1     0     0 ;
    %  1     2     0     0     0 ;
    %  1     0     0     1     0 ;
    %  0     1     1     0     0 ;
    %  0     0     0     0     1 ]
    %>> sum(Output .* [1:5],2)
    % [5 ;
    %  5 ;
    %  5 ;
    %  5 ;
    %  5 ;
    %  5 ;
    %  5 ]
    % Another Example:
    % Base = 2 ;
    % Target = 11 ;
    % Output = NumberDistributionStates(Base,Target)
    % >> Output :
    % [11    0 ;
    %  9     1 ;
    %  7     2 ;
    %  5     3 ;
    %  3     4 ;
    %  1     5 ]
    %>> sum(Output .* [1:2],2)
    % [11 ;
    %  11 ;
    %  11 ;
    %  11 ;
    %  11 ;
    %  11 ]

    if (floor(Base)==ceil(Base)) && (floor(Target)==ceil(Target)) %&& (Base>1) && (Target>1)
        % Initial
        Cellar = {1 , zeros(1,Base)} ; Cellar{1,2}(1,1) = 1 ;
        % Next Alorithm
        for i = 2 : Target
            Casper = Chromosome(i) ;
            CalU = [] ;
            for jj = 1:size(Casper,1)
                Array = Combiner(Cellar{Casper(jj,1),2},Cellar{Casper(jj,2),2}) ;
                CalU = [CalU ; Array] ;
            end
            if i <= Base
                Appi = zeros(1,Base) ; Appi(i) = 1 ;
                CalU = [CalU ; Appi] ;
            end
            CalU = flipud(unique(CalU,'rows')) ;
            Cellar(i,:) = {i , CalU} ;
        end
        if Target == 1
          CalU = Cellar{1,2} ;
        end
    else
        error('Error, Choose two integers greater than 1 numbers.')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Son = Chromosome(I)
    a = I-1 ;
    b = 1 ;
    Son = [] ;
    while a >= b
        Kid = [a b] ;
        Son = [Son ; Kid] ;
        a = a-1 ;
        b = b+1 ;
    end
end

function CO = Combiner(A , B)
    CO = [] ;
    for i = 1:size(A,1)
        for j = 1:size(B,1)
            Cal = A(i,:) + B(j,:) ;
            CO = [CO ; Cal] ;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
