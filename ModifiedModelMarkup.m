%% modified model

% run this section will gives the upstream markup and downstream markup for
% my modified model

% given initial downstream m^d
% calculate passthrough effect dpdu

% update upstream markup and down from DFOC and UFOC

%load data

%load('data632022.mat')
tic
downmarkup_n666=zeros(size(pr));
upmarkup_n666=zeros(size(pr));
for i =1:height(uniqm)

    if mod(i,1000)==0
        disp(i);%track progress
    end
products_m=productsall(ic==i,:);
t_codeshare_m=t_codeshare(ic==i);
nofcodeshare_m=sum(t_codeshare_m);

if nofcodeshare_m>0
    
    ti_m=ti(ic==i);%ticketing carrier 
    op_m=op(ic==i);%operating carrier 
    
    %recode virtual codeshare product due to limited pricing passthrough
    %effect

    op_m(products_m.Codeshare=="1"&products_m.Newonline==1)=...
        ti(products_m.Codeshare=="1"&products_m.Newonline==1);

    op_m_c=op_m(t_codeshare_m==1);%codeshare product's operating carrier
    T=1*(ti_m==ti_m');%ticketing carrier ownership matrix
    T_u=1*(op_m==op_m');%operating carrier ownership matrix
    
    %operating carrier ownership matrix indexed by the ticketing carrier
    T_u_d=1*(ti_m==op_m');

    %operating carrier ownership matrix indexed by the ticketing carrier
    T_d_u=1*(op_m==ti_m');

    Omega_d_m=T.*dsdp_nest{i};%price derivative with ownership matrix, transposed 
    
    pr1000_m=pr1000(ic==i);%price
    s_m=share(ic==i);%share
    J=size(s_m,1);
    within=withinshare(ic==i);
    
    %nestp=0;
    nestp=betan(end);
    delta_nest_m=delta(ic==i)-nestp*log(within);
    Dg=sum(exp(delta_nest_m/(1-nestp)));
    longterm=exp(delta_nest_m/(1-nestp))*(1+(nestp*Dg^(nestp-1)))/(Dg+Dg^nestp);
    dsdd=1/(1-nestp)*(diag(s_m)-s_m*longterm');
    
    %initial Downstream markup from simple upstream downstream model 
    downmarkup_m = -Omega_d_m\s_m;
    diff = 1;
    while diff>1e-10
    
        %compute the second derivative
        term1=s_m*s_m';
        term1remap=repmat(term1,1,1,J);%third dimensional m 
        Dg=sum(exp(delta_nest_m/(1-nestp)));
        
        dsdd2_1=nestp/(1-nestp)*Dg^(nestp-2)*term1remap.*(permute(exp(delta_nest_m/(1-nestp)),[3,2,1]));
        dsdd2_2=-1/(1-nestp)*(nestp*Dg^(nestp-1)+1)*repmat(permute(dsdd,[1,3,2]),1,J,1).*repmat(s_m',J,1,J);
        dsdd2_3=-1/(1-nestp)*(nestp*Dg^(nestp-1)+1)*repmat(permute(dsdd,[3,1,2]),J,1,1).*repmat(s_m,1,J,J);
        
        dsddown_j_j_m=zeros(J,J,J);
        dsddown_j_j_m(bsxfun(@plus,[1:J+1:J*J]',[0:J-1]*J*J))=1/(1-nestp)*dsdd;
        
        dsdp2=betan(1)^2*(dsdd2_1+dsdd2_2+dsdd2_3+dsddown_j_j_m);
        T3=repmat(T,1,1,J); 
        % compute the price pass-through
        
        G=dsdp_nest{i}'+squeeze(sum(T3.*dsdp2.*downmarkup_m,1))' + T.*dsdp_nest{i};
        % create the codeshare matrix I first
        I_c=op_m==op_m_c';II_c=ti_m==op_m_c';
        H=(I_c-II_c).*dsdp_nest{i}(:,t_codeshare_m==1);
        
        % upstream downstream price pass-through
        dudd=G\H; %(j,k) is dpj/dpk

        
        
        % solve for system of equations combining two FOC
        % dsdu
        
        dsdp_u=(dudd'*dsdp_nest{i});
        
        
        FOC11=s_m;
        FOC21=s_m(t_codeshare_m==1)+T_d_u(t_codeshare_m==1,:).*dudd'*s_m;
        
        FOC12=T.*dsdp_nest{i};
        FOC22=T_d_u(t_codeshare_m==1,:).*dsdp_u;


        FOC13=(T_u_d.*dsdp_nest{i});
        FOC13=FOC13(:,t_codeshare_m==1);

        FOC23=T_u(t_codeshare_m==1,:).*dsdp_u;
        FOC23=FOC23(:,t_codeshare_m==1);


        Markup=-[FOC12, FOC13;FOC22,FOC23]\[FOC11;FOC21];
        downmarkup_m_old=downmarkup_m;
        downmarkup_m=Markup(1:end-nofcodeshare_m);
        upmarkup_m=Markup(end-nofcodeshare_m+1:end);
        diff=abs(max(downmarkup_m-downmarkup_m_old));
    end
%wanna see the difference between modified model and Gayle's model
%a=[downmarkup_m,-Omega_d_m\s_m];
upmarkup_n666(ic==i&t_codeshare==1)=upmarkup_m;
downmarkup_n666(ic==i)=downmarkup_m;
else 
    ti_m=ti(ic==i);%ticketing carrier 
    T=1*(ti_m==ti_m');%ticketing carrier ownership matrix
    Omega_d_m=T.*dsdp_nest{i};%price derivative with ownership matrix, transposed 
    s_m=share(ic==i);%share

    downmarkup_n666(ic==i)=-Omega_d_m\s_m;
end 

end

totalmarkup666=downmarkup_n666+upmarkup_n666;
mc_n_666=pr-totalmarkup666*1000;
toc

%%




