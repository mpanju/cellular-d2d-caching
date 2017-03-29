clear all
close all
set(0,'RecursionLimit',1000)

cont_m=10; %%Content for which the P_hit has to be derived
samp=10000;%%NO of samples of X and X_ON
M=100; %%Total Contents
Gam=2;
lam=1;%rate of request
P_pop=(1./([1:M].^Gam))/sum((1./([1:M].^Gam)));
% P_pop=ones(1,M).*1/M;
plot([1:M],P_pop);
CDF_P=P_pop(1)*ones(1,M);
for i=2:M
    CDF_P=CDF_P+[zeros(1,i-1) P_pop(i)*ones(1,M-(i-1))];
end;

% figure
% plot([1:M],CDF_P);
% hist(gen_cdf_dis([1:M],CDF_P,1,10000),[1:M]);


%[sort_time_line, req_order] = requests( lam, P_pop, 1000, cont_m);
% a=load('requestorder.mat');
% b=load('timeline.mat');
% req_order=a.req_order;
% sort_time_line=b.sort_time_line;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_samp=100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sort_time_line=exprnd(1,1,no_samp);
% tmp=[];
%     for j=1:no_samp
%         tmp(j)=sum(sort_time_line(1:j));
%     end
sort_time_line=cumsum(sort_time_line);
req_order=gen_cdf_dis(1:M,CDF_P,1,no_samp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
%plot(sort_time_line, req_order)
% t_end=max(find(req_order==1));
t_end=no_samp;
% C_size=1:5:100; %%for cont 10
% C_size=1:100:1000;%%for cont 100
%C_size=1:20:200;%%for cont 100
C_size=1:5:21;
% C_size=1:2:40;%%for cont 1
P_hit=[];
P_hit_che=[];
Geo_mean=[];
T_C=[];
T_C_rand_vec=[];
Miss_time=[];
var_T_C_rand=[];
mean_T_C_rand=[];
P_hit_TTL=[];
t_phit_lru=[];
t_phit_ttl=[];
t_phit_ttl_ex=[];
for q=C_size
    C=q;
    [req_root, time_root, Miss_time, T_C_rand_vec, sum_vec, sum_mean, Geo_mean_tmp,P_hit_tmp] = LRU_ON_OFF( lam, cont_m, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C );
    

    Timer_c=Che_T_C( M, C, C, P_pop, lam,1e2, 10e6);    
    [TTL_req_root, TTL_time_root, TTL_Miss_time_2nd_level, TTL_T_C_rand_vec_2nd_level, TTL_sum_vec_2nd_level, TTL_sum_mean_2nd_level, TTL_Geo_mean_tmp_2nd_level,TTL_P_hit_tmp_2nd_level] = TTL_ON_OFF( lam, cont_m, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C, Timer_c );
    P_hit_TTL=[P_hit_TTL TTL_P_hit_tmp_2nd_level];
    
    P_hit=[P_hit P_hit_tmp];
    Geo_mean=[Geo_mean Geo_mean_tmp];
    
    che_const=Che_T_C( M, C, cont_m, P_pop, lam,1e2, 10e6);
    T_C=[T_C che_const];
     lam_m=lam*P_pop(cont_m);

    %     Geo_p=exp(-lam_m*T_C(end));
%     Geo_p=1/(1+ceil(Geo_mean(end)));
    Geo_p=exp(-lam_m*C/(1-lam_m));
    P_hit_che=[P_hit_che P_hit_approx(lam,P_pop,cont_m, T_C(end), Geo_p)];
    
    
    
   

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Random sum and T_C %%%%%%%%%%%%%%%%%%%%%%
%     figure;
%     plot(0:10*sum_mean,expcdf(0:10*sum_mean,sum_mean),'color','r')
%     hold
%     plot(0:10*sum_mean,expcdf(0:10*sum_mean,((-1/lam_m)/log(1-Geo_p))),'color','b')
%     
%     [N,X_bin]=hist(sum_vec,max([100,length(sum_vec)/10]));
%     plot(X_bin,cumsum(N/length(sum_vec)),'color','g');
%     
%     
%     
%     figure;
%     hist(T_C_rand_vec,100);
    var_T_C_rand=[var_T_C_rand var(T_C_rand_vec)];
    mean_T_C_rand=[mean_T_C_rand mean(T_C_rand_vec)];
    %%%%%%%%%%%%%%%%%%%%%      Miss Histogram      %%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    [N_miss,X_miss_bin]=hist(Miss_time,max([100,length(Miss_time)/10]));
    
    hist(Miss_time);
    
    
    figure;
    plot(X_miss_bin,cumsum(N_miss/length(Miss_time)),'color','g');
    
    pdf_1_rs_che=exppdf(0:max(Miss_time),((-1/lam_m)/log(1-exp(-che_const*lam_m))));
    pdf_1=exppdf(0:max(Miss_time),((-1/lam_m)/log(1-exp(-C*lam_m/(1-lam_m)))));
    pdf_2=exppdf(0:max(Miss_time),1/lam_m);
    pdf_3=[zeros(1,floor(C/(1-lam_m))) 1 zeros(1,floor(max(Miss_time))-floor(C/(1-lam_m)))];
    hold;
    
    LAM1=((-1/lam_m)/log(1-exp(-che_const*lam_m)));
    LAM2=1/lam_m;
    K1=LAM1*LAM2/(LAM2-LAM1);
    K2=LAM1*LAM2/(LAM1-LAM2);
    
    
    pdf_miss=conv(pdf_1,conv(pdf_2,pdf_3));
    pdf_miss_rs_che=conv(pdf_1_rs_che,conv(pdf_2,pdf_3));
%     T_pdf=che_const:max(Miss_time);
    T_pdf=0:max(Miss_time);
%     pdf_miss=[zeros(1,floor(che_const)) K1*exp(-LAM1*(T_pdf-che_const))+K2*exp(-LAM2*(T_pdf-che_const))];
    plot(0:length(pdf_miss)-1,cumsum(pdf_miss),'color','r');
    plot(0:length(pdf_miss_rs_che)-1,cumsum(pdf_miss_rs_che),'color','c');
    
    
    % Che's Miss rate/2nd level request rate
    
    SIGMA_0 = (1/lam_m)*exp(lam_m*che_const);
    
    che_pdf_form=SIGMA_0*exp(-SIGMA_0.*(T_pdf-che_const));
    che_cdf_form=1-exp(-SIGMA_0.*(T_pdf-che_const));
    che_truncation=[zeros(1,max(find(T_pdf<che_const))) ones(1,length(T_pdf)-max(find(T_pdf<che_const)))];
    che_cdf_form=che_cdf_form.*che_truncation;
%     che_pdf_form=[zeros(1,floor(che_const)) exppdf(((1/lam_m)-che_const), 0:length(T_pdf)-1)];
%     che_pdf_form=exppdf(((1/lam_m)*exp(lam_m*che_const)), 0:length(T_pdf)-1);
%     plot(T_pdf,cumsum(che_pdf_form),'color','b');
    %plot(T_pdf,(che_cdf_form),'color','b');
    plot([0:floor(che_const) (1:10*SIGMA_0)+floor(che_const)],[zeros(1,floor(che_const)) expcdf(0:10*SIGMA_0,SIGMA_0)] ,'color','b');

    total_phit_lru=0;
    total_phit_ttl=0;
    total_phit_ttl_ex=0;
    for cm=1:M

        %[req_root, time_root, Miss_time, T_C_rand_vec, sum_vec, sum_mean, Geo_mean_tmp,P_hit_tmp] = LRU_ON_OFF( lam, cm, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C );
        Timer_c=Che_T_C( M, C, C, P_pop, lam,1e2, 10e6);    
        [TTL_req_root, TTL_time_root, TTL_Miss_time_2nd_level, TTL_T_C_rand_vec_2nd_level, TTL_sum_vec_2nd_level, TTL_sum_mean_2nd_level, TTL_Geo_mean_tmp_2nd_level,TTL_P_hit_tmp_2nd_level] = TTL_ON_OFF( lam, cm, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C, .1*Timer_c );
        %[TTL_req_root, TTL_time_root, TTL_Miss_time_2nd_level, TTL_T_C_rand_vec_2nd_level, TTL_sum_vec_2nd_level, TTL_sum_mean_2nd_level, TTL_Geo_mean_tmp_2nd_level,TTL_EX_P_hit_tmp_2nd_level] = TTL_ON_OFF_exp_phit( lam, cm, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C, Timer_c );
        [TTL_req_root, TTL_time_root, TTL_Miss_time_2nd_level, TTL_T_C_rand_vec_2nd_level, TTL_sum_vec_2nd_level, TTL_sum_mean_2nd_level, TTL_Geo_mean_tmp_2nd_level,TTL_EX_P_hit_tmp_2nd_level] = TTL_ON_OFF( lam, cm, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C, .5*Timer_c );
        [req_root, time_root, Miss_time, T_C_rand_vec, sum_vec, sum_mean, Geo_mean_tmp,P_hit_tmp] = TTL_ON_OFF( lam, cm, P_pop, CDF_P, sort_time_line(1:t_end), req_order(1:t_end), C, Timer_c );
        total_phit_lru=total_phit_lru+P_pop(cm)*P_hit_tmp;
        total_phit_ttl=total_phit_ttl+P_pop(cm)*TTL_P_hit_tmp_2nd_level;
        total_phit_ttl_ex=total_phit_ttl_ex+P_pop(cm)*TTL_EX_P_hit_tmp_2nd_level;
        cm
    end;

    t_phit_lru=[t_phit_lru total_phit_lru];
    t_phit_ttl=[t_phit_ttl total_phit_ttl];
    t_phit_ttl_ex=[t_phit_ttl_ex total_phit_ttl_ex];
q
end
figure;
semilogx(C_size,P_hit,'color','k');

hold
semilogx(C_size,P_hit_TTL,'color','y');
semilogx(C_size,Geo_mean./(1+Geo_mean),'color','r');
P_hit_valentina=1-exp(-lam*P_pop(cont_m).*T_C);
semilogx(C_size,P_hit_valentina,'color','b');
semilogx(C_size,P_hit_che,'color','g');
semilogx(C_size,1./(1+Geo_mean),'color','c');


figure;
semilogx(C_size,mean_T_C_rand,'color','b');
hold;
semilogx(C_size,T_C,'color','g');
semilogx(C_size,C_size./(1-lam_m),'color','c');


semilogx(C_size,sqrt(var_T_C_rand)./mean_T_C_rand,'color','r');




figure;
semilogx(C_size,t_phit_lru,'color','k');

hold
semilogx(C_size,t_phit_ttl,'color','r');
semilogx(C_size,t_phit_ttl_ex,'color','b');