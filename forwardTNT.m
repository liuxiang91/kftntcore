function [ outData ,dataf] = forwardTNT(A,C,Q,R,INITV,INITX,regCoef,b,rho,data ,predictionLength)
% A,C,Q,R,INITV,INITX : kalman filter parameters
% rho: coffidence region size
% b: threshold of test
% regCoef: logistic regression coefficient
% predictionLength: number of periods into the future

dataf={};
rhoOld=rho;
bOld=b;


outData={'id','Visit','MD','IOP','PSD','FilteredMD','FilteredIOP','FilteredPSD','FilteredMD SD','FilteredIOP SD','FilteredPSD SD','Probability of Progression','TNT Test'};


numPat=size(data,1)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
probOfProg=cell(numPat+1,1);
ifTest=cell(numPat+1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate TNT for every single patient;
for i = 2:numPat+1
    rho=rhoOld;
    b=bOld;
    obs_data = data{i,3};
    trueCurrent=size(obs_data,2);
    obs_data=[obs_data nan(9,predictionLength)];
    T = size(obs_data,2);
    
    Baseline_mat = repmat(mean(obs_data(1:3,trueCurrent-1:trueCurrent)')',1,T);
    age= [data{i,5} data{i,5}(trueCurrent)+180:180:(data{i,5}(trueCurrent)+180*(predictionLength)) ];
    race = ones(1,T)*data{i,6}(end);
    sex=ones(1,T)*data{i,7}(end);
    Baseline_MD = Baseline_mat(1,:);
    Baseline_IOP = Baseline_mat(2,:);
    Baseline_PSD = Baseline_mat(3,:);
    probOfProg{i-1} = zeros(T,1);
    ifTest{i-1} = zeros(T,1);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start Prediction
    current_time = trueCurrent;
    fprintf('--start of the patient--\n');
    while current_time < T
        [Xfilt, Vfilt, ~, ~] = kalman_filter(obs_data(:,1:current_time),A,C,Q,R,INITX,INITV);
        V=Vfilt(:,:,current_time);
        t = 1;
        fprintf('%d=Current time, Current Cov=%f\n',current_time,norm(V));
        while t <= (T-current_time)
            current_x = Xfilt(:,current_time);
            predict_x = A^t*current_x;
            accum = 0;
            for j = 0:t-1
                accum = accum + A^j*Q*A^j';
            end
            predict_v = (A^t*V*A^t' + accum);
            try
            a=regCoef(2:10)' ;
            wx = a*predict_x + regCoef(11:13)'*[Baseline_MD(current_time+t) Baseline_IOP(current_time+t) Baseline_PSD(current_time+t)]' +regCoef(1)+...
                regCoef(14:16)'*[age(current_time+t) race(current_time+t) sex(current_time+t)]' ;
            prob=1/(1+exp(-(wx+ sqrt(chi2inv(rho,9)*a*predict_v*a'))));



            fprintf('%d Predict %d step\n',current_time+t,t);
            fprintf('  Cov=%f\n',norm(predict_v));
            fprintf('  wx=%f, conf=%f\n',wx,sqrt(chi2inv(rho,9)*a*predict_v*a'));
            fprintf('  Mean prog prob = %f, with c.i. = %f\n',1/(1+exp(-(wx))),prob);
            catch
            prob=0.5;
            end
            probOfProg{i-1}(current_time+t) = prob;
            %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Call a test
            
            if prob >= b || t>3
                if t>3
                    fprintf('  Test because too long\n');
                end
                if prob>=b
                    b=0.4; %once progressed change to a more aggresive setting
                    rho=0.1;
                    fprintf('  Test because above threshold\n')
                end
                ifTest{i-1}(current_time+t) = 1;
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
    [X, V, ~, ~] = kalman_filter(obs_data,A,C,Q,R,INITX,INITV);
    temp=[];
    for k=1:trueCurrent-1
        outData=[outData; {data{i,1}  data{i,2}(1,k) data{i,3}(1,k) data{i,3}(2,k) data{i,3}(3,k) X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k)) 'Past' 'Past' }];
        temp=[temp;  [  data{i,2}(1,k) data{i,3}(1,k) data{i,3}(2,k) data{i,3}(3,k) X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k)) 0 nan] ];
        
    end
    k=trueCurrent;
    outData=[outData; {data{i,1}  data{i,2}(1,k) data{i,3}(1,k) data{i,3}(2,k) data{i,3}(3,k) X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k))  'Current' 'Current' }];
    temp=[temp; [ data{i,2}(1,k) data{i,3}(1,k) data{i,3}(2,k) data{i,3}(3,k) X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k)) 0 nan ]];
    for k=trueCurrent+1:T
        outData=[outData; {data{i,1}  outData{end,2}+6 nan nan nan X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k))  probOfProg{i-1}(k) ifTest{i-1}(k) }];
        temp=[temp; [  temp(end,1)+6 nan nan nan X(1,k)  X(2,k)  X(3,k) sqrt(V(1,1,k)) sqrt(V(2,2,k)) sqrt(V(3,3,k))  ifTest{i-1}(k) probOfProg{i-1}(k)]];
        
    end
    dataf=[dataf; temp];
end
end

