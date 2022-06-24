//#include <Rcpp.h>
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;
using namespace R;
// [[Rcpp::depends(RcppArmadillo)]]


arma::mat outer_dif(double nt,double nc);
arma::mat outer_sd(double nt,double nc);
arma::vec match_cpp(arma::vec v1, arma::vec v2, double rev);
int find_min(arma::vec a);
List individual_p_cpp(double xt, double xc, double nt, double nc,double pcupper,double pclower,
                     arma::vec theta_grd, double Pgrid,
                     double dev_region_num);

// [[Rcpp::export]]
List ExactGen(arma::vec Xt, arma::vec Xc,
              arma::vec Nt, arma::vec Nc,
              arma::vec pcupper,arma::vec pclower,
              double Thetagrid,
              double Pgrid,
              double cov_prob,
              double B,
              double dev_region_num){
  int K = Xt.n_rows;
  //MH
  vec pkt = Xt/Nt;
  vec pkc = Xc/Nc;
  vec theta_estmh = pkt-pkc;
  vec theta_varmh = (pkt%(1-pkt))/Nt + (pkc%(1-pkc))/Nc;
  vec weightmh = (Nt%Nc/(Nt+Nc))/(sum((Nt%Nc/(Nt+Nc))));
  double MH_org = sum(theta_estmh%weightmh);
  double MH_sd = sqrt(sum(theta_varmh%weightmh%weightmh)); 
  double MH_CIL = MH_org - R::qnorm5((1+cov_prob)/2,0,1,1,0)*MH_sd;
  double MH_CIU = MH_org + R::qnorm5((1+cov_prob)/2,0,1,1,0)*MH_sd;
  double mhpchiq = (MH_org*MH_org)/(MH_sd*MH_sd);
  double MH_p = 1-R::pchisq(mhpchiq,1,1,0);
  
  //MH-0.5 namely co===correction.
  vec Xtco(K,fill::zeros);
  vec Xcco(K,fill::zeros);
  vec Ntco(K,fill::zeros);
  vec Ncco(K,fill::zeros);
  
  for(int i=0;i<K;i++){
    if((pkt(i)==0)|(pkc(i)==0)){
      Xtco(i) = Xt(i)+0.5;
      Xcco(i) = Xc(i)+0.5;
      Ntco(i) = Nt(i)+1;
      Ncco(i) = Nc(i)+1;
    }
    else{
      Xtco(i) = Xt(i);
      Xcco(i) = Xc(i);
      Ntco(i) = Nt(i);
      Ncco(i) = Nc(i);
    }
  }
  
  vec pktco = Xtco/Ntco;
  vec pkcco = Xcco/Ncco;
  vec theta_estmhco = pktco-pkcco;
  vec theta_varmhco = (pktco%(1-pktco))/Ntco + (pkcco%(1-pkcco))/Ncco;
  vec weightmhco = (Ntco%Ncco/(Ntco+Ncco))/(sum((Ntco%Ncco/(Ntco+Ncco))));
  double MH_orgco = sum(theta_estmhco%weightmhco);
  double MH_sdco = sqrt(sum(theta_varmhco%weightmhco%weightmhco)); 
  double MH_CILco = MH_orgco - R::qnorm5((1+cov_prob)/2,0,1,1,0)*MH_sdco;
  double MH_CIUco = MH_orgco + R::qnorm5((1+cov_prob)/2,0,1,1,0)*MH_sdco;
  double mhpchiqco = (MH_orgco*MH_orgco)/(MH_sdco*MH_sdco);
  double MH_pco = 1-R::pchisq(mhpchiqco,1,1,0);
  
  //Exact Method
  vec thetazero_temp={abs(MH_CIL),
                      abs(MH_CIU),
                      abs(MH_CILco),
                      abs(MH_CIUco)};
  double thetazero = arma::max(thetazero_temp);
  vec thetagrad_temp1 = {1,5*thetazero};
  vec thetagrad_temp2 = {1,5*thetazero};
  double theta_grdL = arma::min(thetagrad_temp1);
  double theta_grdU = arma::min(thetagrad_temp2);
  vec theta_grd(Thetagrid,fill::zeros);
  theta_grd(0) = 0;
  for(int i=1;i<Thetagrid;i++){
    theta_grd(i) = -theta_grdL+i*((theta_grdU+theta_grdL)/(Thetagrid-1));
  }
  std::sort(theta_grd.begin(),theta_grd.end());

  
  mat P_value_up(K,Thetagrid,fill::zeros);
  mat P_value_low(K,Thetagrid,fill::zeros);
  
  
  for(int i=0;i<K;i++){
    //double pkcupper = pcupper(i);
    //double pkclower = pclower(i);
    List P_indi = individual_p_cpp(Xt(i),Xc(i),
                                  Nt(i),Nc(i),
                                  pcupper(i),pclower(i),
                                  theta_grd,Pgrid,dev_region_num);
    vec P_indiup = P_indi["p_up"];
    vec P_indilow = P_indi["p_low"];
    P_value_up.row(i) = P_indiup.t();
    P_value_low.row(i) = P_indilow.t();
    
    
    for(int j=0;j<Thetagrid;j++){
      rowvec max_temp_low = P_value_low(i,span(0,(Thetagrid-j-1)));
      rowvec max_temp_up = P_value_up(i,span(j,Thetagrid-1));
      P_value_low(i,(Thetagrid-j-1))= arma::max(max_temp_low.t());
      P_value_up(i,j) = arma::max(max_temp_up.t());
    }
  }
  //DEBUG
  
  P_value_low = P_value_low/(1+1e-3);
  P_value_up = P_value_up/(1+1e-3);

  //permutation
  //mat tnull(B,3);
  mat permut(K,B,fill::randu);
  permut = permut/(1+1e-3);
  //extend weight (K*1) to (K*B)
  mat weight_extend(K,B,fill::ones);
  for(int i=0;i<B;i++){
    weight_extend.col(i) = weightmh;
  }
  rowvec iden_exact = sum(permut%weight_extend,0);
  rowvec asin_exact = sum(asin(sqrt(permut))%weight_extend,0);

  mat qnorm_exact(K,B,fill::zeros);
  for(int i=0;i<K;i++){
    for(int j=0;j<B;j++){
      qnorm_exact(i,j) = R::qnorm5(permut(i,j),0,1,1,0);
    }
  }
  rowvec invnorm_exact = sum(qnorm_exact%weight_extend,0);


  double alpha0 = (1+cov_prob)/2;
  vec cutoff(3);
  vec P = {0.5,1-alpha0};
  vec Q1 = quantile(iden_exact.t(),P);
  cutoff(0) = Q1(1);
  vec Q2 = quantile(asin_exact.t(),P);
  cutoff(1) = Q2(1);
  vec Q3 = quantile(invnorm_exact.t(),P);
  cutoff(2) = Q3(1);


  //mat t1(Thetagrid,3);
  //mat t2(Thetagrid,3);
  mat weight_extend_exact(K,Thetagrid,fill::ones);
  for(int i=0;i<Thetagrid;i++){
    weight_extend_exact.col(i) = weightmh;
  }
  //p value lower
  rowvec iden_exact_lower = sum(P_value_low%weight_extend_exact,0);
  rowvec asin_exact_lower = sum(asin(sqrt(P_value_low))%weight_extend_exact,0);

  mat qnorm_exact_lower(K,Thetagrid,fill::zeros);
  for(int i=0;i<K;i++){
    for(int j=0;j<Thetagrid;j++){
      qnorm_exact_lower(i,j) = R::qnorm5(P_value_low(i,j),0,1,1,0);
    }
  }
  rowvec invnorm_exact_lower = sum(qnorm_exact_lower%weight_extend_exact,0);
  //p value upper
  rowvec iden_exact_upper = sum(P_value_up%weight_extend_exact,0);
  rowvec asin_exact_upper = sum(asin(sqrt(P_value_up))%weight_extend_exact,0);

  mat qnorm_exact_upper(K,Thetagrid,fill::zeros);
  for(int i=0;i<K;i++){
    for(int j=0;j<Thetagrid;j++){
      qnorm_exact_upper(i,j) = R::qnorm5(P_value_up(i,j),0,1,1,0);
    }
  }
  rowvec invnorm_exact_upper = sum(qnorm_exact_upper%weight_extend_exact,0);


  // 
  // //find exact's CI
  // //identity
  uvec iden_exact_lower_ind = find(iden_exact_lower>cutoff(0));
  double iden_exact_lower_ind_len = iden_exact_lower_ind.n_elem;
  uvec iden_exact_upper_ind = find(iden_exact_upper>cutoff(0));
  double iden_exact_upper_ind_len = iden_exact_upper_ind.n_elem;

  vec iden_theta_grad_min(iden_exact_lower_ind_len);
  vec iden_theta_grad_max(iden_exact_upper_ind_len);
  for(int i=0;i<iden_exact_lower_ind_len;i++){
    iden_theta_grad_min(i) = theta_grd(iden_exact_lower_ind(i));
  }
  for(int i=0;i<iden_exact_upper_ind_len;i++){
    iden_theta_grad_max(i) = theta_grd(iden_exact_upper_ind(i));
  }

  double Iden_CIL = arma::min(iden_theta_grad_min);
  double Iden_CIU = arma::max(iden_theta_grad_max);
  

  // //asin
  uvec asin_exact_lower_ind = find(asin_exact_lower>cutoff(1));
  double asin_exact_lower_ind_len = asin_exact_lower_ind.n_elem;
  uvec asin_exact_upper_ind = find(asin_exact_upper>cutoff(1));
  double asin_exact_upper_ind_len = asin_exact_upper_ind.n_elem;

  vec asin_theta_grad_min(asin_exact_lower_ind_len);
  vec asin_theta_grad_max(asin_exact_upper_ind_len);
  for(int i=0;i<asin_exact_lower_ind_len;i++){
    asin_theta_grad_min(i) = theta_grd(asin_exact_lower_ind(i));
  }
  for(int i=0;i<asin_exact_upper_ind_len;i++){
    asin_theta_grad_max(i) = theta_grd(asin_exact_upper_ind(i));
  }

  double Asin_CIL = arma::min(asin_theta_grad_min);
  double Asin_CIU = arma::max(asin_theta_grad_max);

  //qnorm
  uvec invnorm_exact_lower_ind = find(invnorm_exact_lower>cutoff(2));
  double invnorm_exact_lower_ind_len = invnorm_exact_lower_ind.n_elem;
  uvec invnorm_exact_upper_ind = find(invnorm_exact_upper>cutoff(2));
  double invnorm_exact_upper_ind_len = invnorm_exact_upper_ind.n_elem;

  vec invnorm_theta_grad_min(invnorm_exact_lower_ind_len);
  vec invnorm_theta_grad_max(invnorm_exact_upper_ind_len);
  for(int i=0;i<invnorm_exact_lower_ind_len;i++){
    invnorm_theta_grad_min(i) = theta_grd(invnorm_exact_lower_ind(i));
  }
  for(int i=0;i<invnorm_exact_upper_ind_len;i++){
    invnorm_theta_grad_max(i) = theta_grd(invnorm_exact_upper_ind(i));
  }

  double invnorm_CIL = arma::min(invnorm_theta_grad_min);
  double invnorm_CIU = arma::max(invnorm_theta_grad_max);
  

  // //exact
  // //identity
  vec iden_esti_temp = abs(iden_exact_upper.t()-iden_exact_lower.t());
  int iden_esti_ind = find_min(iden_esti_temp);
  double iden_esti = theta_grd(iden_esti_ind);

  //asin
  vec asin_esti_temp = abs(asin_exact_upper.t()-asin_exact_lower.t());
  int asin_esti_ind = find_min(asin_esti_temp);
  double asin_esti = theta_grd(asin_esti_ind);

  //qnorm
  vec invnorm_esti_temp = abs(invnorm_exact_upper.t()-invnorm_exact_lower.t());
  int invnorm_esti_ind = find_min(invnorm_esti_temp);
  double invnorm_esti = theta_grd(invnorm_esti_ind);
  
  // 
  // //calculate p values
  int n0=0;
  for(int i=0;i<Thetagrid;i++){
    if(theta_grd(i)==0){n0=i;}
    else{continue;}
  }
  //uvec n0=find(theta_grd==0);
  //iden
  double c1_iden_lower = iden_exact_lower(n0);
  double c2_iden_upper = iden_exact_upper(n0);
  uvec c1_iden_temp = iden_exact.t()<=c1_iden_lower;
  uvec c2_iden_temp = iden_exact.t()<=c2_iden_upper;
  vec iden_pvalue_temp1={arma::accu(c1_iden_temp)/B,arma::accu(c2_iden_temp)/B};
  vec iden_pvalue_temp = {1.0,2*arma::min(iden_pvalue_temp1)};
  double iden_pvalue = arma::min(iden_pvalue_temp);

  //asin
  double c1_asin_lower = asin_exact_lower(n0);
  double c2_asin_upper = asin_exact_upper(n0);
  uvec c1_asin_temp = asin_exact.t()<=c1_asin_lower;
  uvec c2_asin_temp = asin_exact.t()<=c2_asin_upper;
  vec asin_pvalue_temp1={arma::accu(c1_asin_temp)/B,arma::accu(c2_asin_temp)/B};
  vec asin_pvalue_temp = {1.0,2*arma::min(asin_pvalue_temp1)};
  double asin_pvalue = arma::min(asin_pvalue_temp);
  //qnorm

  double c1_invnorm_lower = invnorm_exact_lower(n0);
  double c2_invnorm_upper = invnorm_exact_upper(n0);
  uvec c1_invnorm_temp = invnorm_exact.t()<=c1_invnorm_lower;
  uvec c2_invnorm_temp = invnorm_exact.t()<=c2_invnorm_upper;
  vec invnorm_pvalue_temp1={arma::accu(c1_invnorm_temp)/B,arma::accu(c2_invnorm_temp)/B};
  vec invnorm_pvalue_temp = {1.0,2*arma::min(invnorm_pvalue_temp1)};
  double invnorm_pvalue = arma::min(invnorm_pvalue_temp);
  


  //RD---->OR
  //iden
  double iden_or_numer = accu(((1-pkc)%(pkc+iden_esti)%(Nt%Nc))/(Nt+Nc));
  double iden_or_denom = accu(pkc%(1-pkc-iden_esti)%(Nt%Nc)/(Nt+Nc));
  double iden_or = iden_or_numer/iden_or_denom;

  //asin
  double asin_or_numer = accu(((1-pkc)%(pkc+asin_esti)%(Nt%Nc))/(Nt+Nc));
  double asin_or_denom = accu(pkc%(1-pkc-asin_esti)%(Nt%Nc)/(Nt+Nc));
  double asin_or = asin_or_numer/asin_or_denom;
  //qnorm
  double invnorm_or_numer = accu(((1-pkc)%(pkc+invnorm_esti)%(Nt%Nc))/(Nt+Nc));
  double invnorm_or_denom = accu(pkc%(1-pkc-invnorm_esti)%(Nt%Nc)/(Nt+Nc));
  double invnorm_or = invnorm_or_numer/invnorm_or_denom;
  //MH
  double MH_or_numer = accu(((1-pkc)%(pkc+MH_org)%(Nt%Nc))/(Nt+Nc));
  double MH_or_denom = accu(pkc%(1-pkc-MH_org)%(Nt%Nc)/(Nt+Nc));
  double MH_or = MH_or_numer/MH_or_denom;
  //MH-0.5
  double MHco_or_numer = accu(((1-pkc)%(pkc+MH_orgco)%(Nt%Nc))/(Nt+Nc));
  double MHco_or_denom = accu(pkc%(1-pkc-MH_orgco)%(Nt%Nc)/(Nt+Nc));
  double MHco_or = MHco_or_numer/MHco_or_denom;

  //output
  vec est={iden_esti,asin_esti,invnorm_esti,MH_org,MH_orgco};
  vec ppp={iden_pvalue,asin_pvalue,invnorm_pvalue,MH_p,MH_pco};
  vec CIL={Iden_CIL,Asin_CIL,invnorm_CIL,MH_CIL,MH_CILco};
  vec CIU={Iden_CIU,Asin_CIU,invnorm_CIU,MH_CIU,MH_CIUco};
  vec OR={iden_or,asin_or,invnorm_or,MH_or,MHco_or};

  return List::create(Rcpp::Named("Estimate")=est,
                      Rcpp::Named("pvalue")=ppp,
                      Rcpp::Named("CILower")=CIL,
                      Rcpp::Named("CIUpper")=CIU,
                      Rcpp::Named("OR")=OR
  );
}

// [[Rcpp::export]]
List individual_p_cpp(double xt, double xc, double nt, double nc,
                     double pcupper,double pclower, arma::vec theta_grd, double Pgrid,
                     double dev_region_num){
  
  double error = 1e-6;
  double pcobs = xc/nc;
  double ptobs = xt/nt;
  double pcobs_correct = (xc+0.5)/(nc+1.0);
  double ptobs_correct = (xt+0.5)/(nt+1.0);
  //pc_grd: column vector
  arma::vec pc_grd(Pgrid,fill::zeros);
  double pcgrad_portion = (pcupper-pclower)/(Pgrid-1.0);
  vec pvalue1(Pgrid,fill::zeros);
  vec pvalue2(Pgrid,fill::zeros);
  for(double i=0; i<Pgrid; i++){
    pc_grd(i) = pclower + pcgrad_portion*i;
  }
  
  arma::mat Munull = outer_dif(nt,nc);
  arma::mat Sdnull = outer_sd(nt,nc);
  double Muobs = ptobs-pcobs;
  double Sdobs = sqrt((pcobs_correct*(1-pcobs_correct))/nc+(ptobs_correct*(1.0-ptobs_correct))/nt);
  int thetagrad_length = theta_grd.n_rows;
  vec p_low(thetagrad_length,fill::zeros);
  vec p_up(thetagrad_length,fill::zeros);

  for(int j=0;j<thetagrad_length;j++){
    double Tobs = (Muobs-theta_grd(j))/Sdobs;
    
    if((arma::min(pc_grd)+theta_grd(j))>=1.0){
      p_low(j) = 1.0;
      p_up(j) =0.0;
    }
    else if((arma::max(pc_grd)+theta_grd(j))<=0.0){
      p_low(j) =0.0;
      p_up(j) =1.0;
    }
    else{
      for(int b=0;b<Pgrid;b++){
        if(((pc_grd(b)+theta_grd(j))>=0) & ((pc_grd(b)+theta_grd(j))<=1)){
          arma::vec pcb(nc+1,fill::zeros);
          arma::vec ptb(nt+1,fill::zeros);
          double ib = ceil((nc+1.0)*pc_grd(b));
          double jb = ceil((nt+1.0)*(pc_grd(b)+theta_grd(j)));
          
          for(int i1=0;i1<=ib;i1++){
            pcb(i1) = R::dbinom(i1,nc,pc_grd(b),false);
          }
          for(int j1=0;j1<=jb;j1++){
            ptb(j1) = R::dbinom(j1,nt,(pc_grd(b)+theta_grd(j)),false);
          }
          
          int nc_remain = nc-ib;
          int nt_remain = nt-jb;
          
          ib = ib + 1;
          jb = jb + 1;
          
          double sum_devregion=0;
          for(int i=0;i<dev_region_num;i++){
            sum_devregion = sum_devregion + i +1.0;
          }
          
          
          double reg = 1;
          int bk = 0;
          int nc_temp;
          while((reg <= dev_region_num) & (bk == 0)){
            nc_temp = ceil((nc_remain*reg)/sum_devregion);
            if((ib+nc_temp-1)>nc){
              nc_temp = nc-ib+1.0;
              bk = 1;}
            for(int i2=ib;i2<ib+nc_temp;i2++){
              pcb(i2) = R::dbinom(i2,nc,pc_grd(b),false);
            }
            ib = ib + nc_temp;
            reg = reg + 1.0;
            if(pcb(ib-1)<1e-15){break;}
          }
          
          double reg2 = 1;
          int bk2 = 0;
          int nt_temp;
          while((reg2 <= dev_region_num) & (bk2 == 0)){
            nt_temp = ceil((nt_remain*reg2)/sum_devregion);
            if((jb+nt_temp-1)>nt){
              nt_temp = nt-jb+1.0;
              bk2 = 1;}
            for(int j2=jb;j2<jb+nt_temp;j2++){
              ptb(j2) = R::dbinom(j2,nt,(pc_grd(b)+theta_grd(j)),false);
            }
            jb = jb + nt_temp;
            reg2 = reg2 + 1.0;
            if(ptb(jb-1)<1e-15){break;}
          }
          vec pc_main(ib,fill::zeros);
          vec pt_main(jb,fill::zeros);
          
          for(int i3=0;i3<ib;i3++){
            pc_main(i3) = pcb(i3);
          }
          for(int j3=0;j3<jb;j3++){
            pt_main(j3) = ptb(j3);
          }

          std::sort(pc_main.begin(),pc_main.end());
          std::sort(pt_main.begin(),pt_main.end());
          
          vec pc_grdcum = cumsum(pc_main);
          arma::uvec pc_grdcum_ind = find(pc_grdcum >= error);
          double pc_grdcum_ind_l = pc_grdcum_ind.n_rows;
          vec pc_main1(pc_grdcum_ind_l,fill::zeros);
          for(int i4=0;i4<pc_grdcum_ind_l;i4++){
            pc_main1(i4) = pc_main(pc_grdcum_ind(i4));
          }
          
          vec pt_grdcum = cumsum(pt_main);
          arma::uvec pt_grdcum_ind = find(pt_grdcum >= error);
          double pt_grdcum_ind_l = pt_grdcum_ind.n_rows;
          vec pt_main1(pt_grdcum_ind_l,fill::zeros);
          for(int i5=0;i5<pt_grdcum_ind_l;i5++){
            
            pt_main1(i5) = pt_main(pt_grdcum_ind(i5));
          }
          
          
          
          //match function
          vec idc = match_cpp(pc_main1,pcb,1);
          vec idt = match_cpp(pt_main1,ptb,1);
          
          int idc_l = idc.n_rows;
          int idt_l = idt.n_rows;
          vec pc_joint(idc_l,fill::zeros);
          vec pt_joint(idt_l,fill::zeros);
          
          for(int i6=0;i6<idc_l;i6++){
            pc_joint(i6) = pcb(idc(i6));
          }
          for(int i7=0;i7<idt_l;i7++){
            pt_joint(i7) = ptb(idt(i7));
          }
          
          mat Joint = pc_joint*pt_joint.t();
          
          
          
          arma::mat Tnull(idc_l,idt_l,fill::zeros);
          for(int i8=0;i8<idc_l;i8++){
            for(int i9=0;i9<idt_l;i9++){
              Tnull(i8,i9) = (Munull(idc(i8),idt(i9))-theta_grd(j))/Sdnull(idc(i8),idt(i9));
            }
          }
          
          
          umat part_gr_mat = Tnull>Tobs;
          mat part_gr_mat1 = Joint%part_gr_mat;
          double part_gr_sum = accu(part_gr_mat1);
          umat part_eq_mat = Tnull==Tobs;

          mat part_eq_mat1 = Joint%part_eq_mat;
          double part_eq_sum = accu(part_eq_mat1);
          double part_lw_sum = accu(Joint)-part_eq_sum-part_gr_sum;

          
          pvalue1(b) = part_gr_sum + 0.5*part_eq_sum;
          pvalue2(b) = part_lw_sum + 0.5*part_eq_sum;
          
        }
        p_up(j) = arma::max(pvalue2);
        p_low(j) = arma::max(pvalue1);
      }
    }
    
    
    
    
  }
  return List::create(Rcpp::Named("p_up")=p_up,
                      Rcpp::Named("p_low")=p_low
  );
}


arma::vec match_cpp(vec v1, vec v2, double rev){
  double l1 = v1.n_rows;
  double l2 = v2.n_rows;
  vec idx(l1,fill::zeros);
  for(int i=0;i<l1;i++){
    for(int j=0;j<l2;j++){
      if(v1(i)==v2(j)){idx(i)=j;}
      else{continue;}
    }
  }
  if(rev==0){return idx;}
  else{return arma::reverse(idx);}
}

arma::mat outer_dif(double nt,double nc){
  arma::mat A(nc+1,nt+1,fill::zeros);
  for(int i=0; i<=nc; i++){
    for(int j=0; j<=nt ; j++){
      A(i,j) = j/nt-i/nc;
    }
  }
  return A;
}

arma::mat outer_sd(double nt,double nc){
  arma::mat A(nc+1,nt+1,fill::zeros);
  int i;
  int j;
  for(i=0; i<=nc; i++){
    for(j=0; j<=nt ; j++){
      double pkt = (j+0.5)/(nt+1.0);
      double pkc = (i+0.5)/(nc+1.0);
      A(i,j) = sqrt((pkt)*(1.0-pkt)/nt+(pkc)*(1.0-pkc)/nc);
    }
  }
  return A;
}

//credit to gallery.rcpp.org/articles/vector-minimum/

int find_min(arma::vec a){
  NumericVector::iterator it = std::min_element(a.begin(),a.end());
  return it-a.begin();
}
