% RECOMPUTE THE MARGIN
% Draft in 5.24.2022
i=4;%test market index 4
test=productsall(ic==i,:); % test market sample
%nestp=0; % test nest parameter
    t_codeshare_m=t_codeshare(ic==i);
    online_m=online(ic==i);
    ti_m=ti(ic==i);%ticketing carrier   
    ti_nc=ti_m(t_codeshare_m==1);
    T=1*(ti_m==ti_m');%ownership matrix
    T_nc=1*(ti_nc==ti_nc');
    index_have_margin=(t_codeshare_m==1);
    Omega_d_m=T.*dsdp_nest{i};%price derivative with ownership matrix, transposed 
    dsdp_nc=dsdp_nest{i}(t_codeshare_m==1,t_codeshare_m==1);
    pr1000_m=pr1000(ic==i);%price
    s_m=share(ic==i);%share
    J=size(s_m,1);
    within=withinshare(ic==i);
    
    s_nc=s_m(t_codeshare_m==1);
    pr1000_nc=pr1000_m(index_have_margin);
    %nestp=0;
    nestp=betan(end);
    delta_nest_m=delta(ic==i)-nestp*log(within);
    Dg=sum(exp(delta_nest_m/(1-nestp)));
    
    longterm=exp(delta_nest_m/(1-nestp))*(1+(nestp*Dg^(nestp-1)))/(Dg+Dg^nestp);
    dsdd=1/(1-nestp)*(diag(s_m)-s_m*longterm');
    %dsdd=dsdp_nest{i}/betan(1);
    
    term1=s_m*s_m';
    term1remap=repmat(term1,1,1,J);%third dimensional m 
    Dg=sum(exp(delta_nest_m/(1-nestp)));
    
    dsdd2_1=nestp/(1-nestp)*Dg^(nestp-2)*term1remap.*(permute(exp(delta_nest_m/(1-nestp)),[3,2,1]));
    dsdd2_2=-1/(1-nestp)*(nestp*Dg^(nestp-1)+1)*repmat(permute(dsdd,[1,3,2]),1,J,1).*repmat(s_m',J,1,J);
    dsdd2_3=-1/(1-nestp)*(nestp*Dg^(nestp-1)+1)*repmat(permute(dsdd,[3,1,2]),J,1,1).*repmat(s_m,1,J,J);
    
    dsddown_j_j_m=zeros(J,J,J);
    dsddown_j_j_m(bsxfun(@plus,[1:J+1:J*J]',[0:J-1]*J*J))=1/(1-nestp)*dsdd;
    
    dsdp2=betan(1)^2*(dsdd2_1+dsdd2_2+dsdd2_3+dsddown_j_j_m); %second derivative of market share with respect to price

%{
%wrong previous formula
    %rearrage share vector to J*J*J matrix with each value in 3rd dimension
    sss=repmat(permute(s_m,[3,2,1]),J,J,1);
    SmSjSk=(s_m*s_m').*sss;
    before_smsjsk=(2*(nestp*Dg^(nestp-1)+1)^2/(1-nestp)^2-nestp/(nestp-1)*Dg^(nestp-1)*(Dg^(nestp-1)+1));
    Sj=zeros(J,J,J);
    %Sj([1:J*J+J+1:J*J*J])=s_m-3*s_m.^2;
    %-3*s_m.^2 is already included in following 3 SmSmk
    Sj([1:J*J+J+1:J*J*J])=s_m;
    SmSmSk=zeros(J,J,J);
    SmSmSk(bsxfun(@plus,[J*J-J+1:-J+1:J]',[0:J-1]*J*J))=(s_m*s_m');
    SmSkSm=zeros(J,J,J);
    SmSkSm(bsxfun(@plus,[1:J:J*J]',[0:J-1]*J*J))=(s_m*s_m');
    SmSkSk=zeros(J,J,J);
    SmSkSk(bsxfun(@plus,[1:J]',[0:J-1]*J*J))=(s_m*s_m');
    dsdp2=(betan(1))^2*((1/(1-nestp)^2)*Sj+before_smsjsk*SmSjSk...
        -(1/(1-nestp)^2)*(nestp*Dg^(nestp-1)+1)*(SmSmSk+SmSkSm+SmSkSk));%second derivative to price
%}

    T3=repmat(T,1,1,J);   
    %Downstream markup
    downmarkup_m=-Omega_d_m\s_m;
    %downmarkup_nc=downmarkup_m(index_have_margin);
    downmarkup_n(ic==i)=downmarkup_m;
    G=dsdp_nest{i}'+squeeze(sum(T3.*dsdp2.*downmarkup_m,1)) + T.*dsdp_nest{i} ;
    H=T.*dsdp_nest{i};
    dddu_test=G\H;
    dddu=G\H;  %price pass-through effect du: upstream, dd downstream
    %upstream sale to price derivative dsdu = dddu'*dsdp

    
    op_m=op(ic==i);%operating carrier  
    op_nc=op_m(index_have_margin);
    
    T_op=1*(op_m==op_m');%ownership matrix   
    T_op_nc=1*(op_nc==op_nc');
    dsdp_u=(dudd'*dsdp_nest{i});%the first derivative of market share with respect to upstream pricing after accounting for price pass-through
    Omega_u_m=T_op_nc.*dsdp_u(index_have_margin,index_have_margin);%upstream firm's ownership matrix
    
    %Upstream markup   
    upmarkup_m=-Omega_u_m\s_nc;
    upmarkup_n3(ic==i&t_codeshare==1&integrated==1)=upmarkup_m; 
    
    % This gives upstream markup given the cases that the set of upstream firm and down stream firms are totally two difrent sets or F_u intersect F_d = None. 
    % Next step is recomputing the markup given horizontal & vetical competition at the same time. 
    
    
    





    
