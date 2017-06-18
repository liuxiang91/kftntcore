%Goals:
%   Calculate number of tests per patient, accuracy, diagnostic delay and
%   cost per patient per period for given b, rho, and cost ratio
%Input:
%   Data sheets
function [accuracy,dd,numTestPerYear,testSpacings] = ...
    TNT(A,C,Q,R,INITV,INITX,regCoef,b,rho,data,warmup,cap)


numPat=size(data,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
probprint=cell(numPat+1,1);
predprint=cell(numPat+1,1);
numTest=zeros(numPat,1);
hit=0;
dd=[];
testSpacings=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate TNT for every single patient;
for i = 2:numPat+1
    obs_data = data{i,3};
    T = size(obs_data,2);
    
    Baseline_mat = data{i,4};
    age= data{i,5};
    race = data{i,6};
    sex=data{i,7};
    Baseline_MD = Baseline_mat(1,:);
    Baseline_IOP = Baseline_mat(2,:);
    Baseline_PSD = Baseline_mat(3,:);
    probprint{i-1} = zeros(T,1);
    predprint{i-1} = zeros(T,1);
    predprint{i-1}(1,1)=1;
    if size(data,2)>8
        prog_vect = data{i,12};
    else
        prog_vect = data{i,8};
    end
    if prog_vect(1) == 1
        hit=hit+1;  %tested when patient progressed
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start Prediction
    current_time = warmup;
    if warmup>=T
        continue;
    end
    %fprintf('--start of the patient--\n');
    while current_time < T
        
        [Xfilt, Vfilt, ~, ~] = kalman_filter(obs_data(:,1:current_time),A,C,Q,R,INITX,INITV);
        V=Vfilt(:,:,current_time);
        t = 1;
        %fprintf('%d=Current time\n',current_time);
        while t <= (T-current_time)
            current_x = Xfilt(:,current_time);
            predict_x = A^t*current_x;
            accum = 0;
            for j = 0:t-1
                accum = accum + A^j*Q*A^j';
            end
            predict_v = (A^t*V*A^t' + accum);
            a=regCoef(2:10)' ;
            % make change here to include additional variables
            if ~isempty(regCoef<=19) % base case no japan variables
                wx = a*predict_x + regCoef(11:13)'*[Baseline_MD(current_time+t) Baseline_IOP(current_time+t) Baseline_PSD(current_time+t)]' +regCoef(1)+...
                    regCoef(14:16)'*[age(current_time+t) race(current_time+t) sex(current_time+t)]' ;
            else
                SE= data{i,8};
                CCT = data{i,9};
                AX=data{i,10};
                DH=data{i,11};
                wx = a*predict_x + regCoef(11:13)'*[Baseline_MD(current_time+t) Baseline_IOP(current_time+t) Baseline_PSD(current_time+t)]' +regCoef(1)+...
                    regCoef(14:20)'*[age(current_time+t) race(current_time+t) sex(current_time+t) SE(current_time+t) CCT(current_time+t) AX(current_time+t) DH(current_time+t)]' ;
            end
            prob=1/(1+exp(-(wx+ sqrt(chi2inv(rho,9)*(a*predict_v*a')))));
            %fprintf('%d Predict %d step\n',current_time+t,t);
            %fprintf('  Mean prog prob = %f, with c.i. = %f\n',1/(1+exp(-(wx))),prob);
            probprint{i-1}(current_time+t) = prob;
            %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Call a test
            if cap==0
                c=Inf;
            else
                c=cap;
            end
            if prob >= b || t>c
                if t>3
                    %fprintf('  Test because too long\n');
                end
                if prob>=b
                    %fprintf('  Test because above threshold\n')
                    %                     Baseline_MD(current_time+t+1:end)=;
                end
                testSpacings=[testSpacings t];
                numTest(i-1,1) = numTest(i-1,1) +1;
                predprint{i-1}(current_time+t) = 1;
                sumprog=sum(prog_vect(current_time:current_time+t));
                if sumprog>=1
                    firstprog=find([zeros(1,current_time) prog_vect(current_time+1:current_time+t)],1,'first'); %record the index of the first progression--could be testing at 2nd progression instance
                    %fprintf('  Catch progression at %d\n',firstprog);
                    dd=[dd current_time+t-firstprog]; %count number of periods since first progression until this test
                    if prog_vect(current_time+t) == 1
                        hit=hit+1;  %tested when patient progressed
                    end
                else
                    %fprintf('  No progression so far\n');
                end
                break;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %No test called
            else
                obs_data(:,current_time+t) = NaN; %you didn't observe these values
            end
            t = t + 1;
        end
        current_time = current_time + t;
    end
    %fprintf('--end of the patient--\n');
    numTest(i-1,1)=2*numTest(i-1,1)/(T-1);
    if sum( predprint{i-1})==1
        testSpacings=[testSpacings 4];
    end
end
if size(data,2)>8
    accuracy=hit/sum([(data{2:end,12})]);
    if (sum([(data{2:end,12})]))==0
        %fprintf('No labeled progression in the input file\n');
    end
else
    accuracy=hit/sum([(data{2:end,8})]);
    if (sum([(data{2:end,8})]))==0
        %fprintf('No labeled progression in the input file\n');
    end
end

numTestPerYear=mean(numTest);
dd=mean(dd);



