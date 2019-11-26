
void StartStateEvolution(function2D<NState>& state_evolution_left,
			 function2D<NState>& state_evolution_right,
			 const function1D<NState>& praStates,
			 NOperators& Operators, const ClusterData& cluster,
			 function1D<Number>& Trace, function2D<double>& Prob)
{
  if (!Operators.empty()){cerr<<"Operators are not empty!"<<endl;}
  function2D<Number> Projection(Prob.size_N(),common::max_size);
  NState nstate(common::max_size,common::max_size), mstate(common::max_size,common::max_size);
  Number ms=0;
  for (int i=0; i<cluster.size(); i++){
    if (praStates[i].empty()) {
      state_evolution_left[i].resize(0);
      state_evolution_right[i].resize(0);
      continue;
    }
    state_evolution_left[i].resize(1);
    state_evolution_right[i].resize(1);
    
    nstate = praStates[i];
    mstate.Evolve(nstate, cluster,  Operators.exp_(0));
    state_evolution_left(i,0).SetEqual(mstate);
    
    mstate.Evolve(nstate, cluster, Operators.exp_(0));
    state_evolution_right(i,0).SetEqual(mstate);
    
    Trace[i] = state_evolution_left[i][0].Project_to_Pra(praStates[i], Projection[i]);
    ms += Trace[i];
  }
  for (int i=0; i<cluster.size(); i++)
    for (int l=0; l<Prob[i].size(); l++)
      Prob[i][l] = divide(Projection[i][l],ms);
}




void Get_N_at_Operator2(function1D<double>& Njc, int ops, const function1D<Number>& Trace, const Number& ms,
			const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right,
			const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{// This routine returns the value of operator N measured at an existing creating operator.
 // If we have diagram with the following configuration <...... psi^+(t_s) ........>
 // We evaluate the following quantity:   <...... psi^+(t_s)psi^+{a}(t_s)psi{b}(t_s) ........>/<...... psi^+(t_s) ........>
 // for two baths {a} and {b}, which have a common index ii.
  Njc=0.0;
  static NState lstate(common::max_size,common::max_size), rstate(common::max_size,common::max_size), lpstate(common::max_size,common::max_size);
  static Number mm;
  for (int ist=0; ist<praStates.size(); ist++){// op_i is the earliest operator inserted - from this point on, the evolution needs to be changed
    double Zi = (Trace[ist]/ms).dbl();
    if (state_evolution_left[ist].size()>Operators.size() && state_evolution_right[ist].size()>Operators.size() && fabs(Zi)>1e-10){
      if (ops!=0) lstate = state_evolution_left[ist][ops-1];
      else lstate = praStates[ist];
      rstate = state_evolution_right[ist][Operators.size()-1-ops];
      int istate = lstate.istate;
      int irstate = rstate.istate;
      if (istate!=irstate) cerr<<"ERROR: states are different"<<istate<<" "<<irstate<<" "<<istep<<endl;
      lpstate.Evolve(lstate, cluster, Operators.exp_(ops));
      //Number mm = lpstate.ScalarProduct(rstate);
      //double Tr_ratio = (mm/Trace[ist]).dbl();
      //if (fabs((Tr_ratio-1.)*Zi)>1e-6) cerr<<"ERROR: Traces["<<ist<<","<<istep<<"]="<<Tr_ratio<<" Zi["<<ist<<"]="<<Zi<<endl;
      for (int jj=0; jj<cluster.Njjs[istate].size(); jj++){
	int ii = cluster.Njjs[istate][jj].ii;
        Multiply(lstate.M, cluster.Njjs[istate][jj].M, lpstate.M);
	lstate.exponent = lpstate.exponent;
	Njc[ii] += (rstate.ScalarProduct(lstate)/ms).dbl();
      }     
    }
  }
}


Number UpdateStateEvolution(function2D<NState>& state_evolution_left, function2D<NState>& state_evolution_right,
			    int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster,
			    function1D<Number>& Trace, function2D<double>& Prob, long long istep, bool cmp_P=true)
{
  static function2D<Number> Projection(Prob.size_N(),common::max_size);
  static NState nstate(common::max_size,common::max_size), mstate(common::max_size,common::max_size);
  Number ms=0;
  Projection=0;

  for (int ist=0; ist<praStates.size(); ist++){// op_i is the earliest operator inserted - from this point on, the evolution needs to be changed
    Trace[ist]=0;
    if (op_i>0 && (state_evolution_left[ist].size()<op_i || state_evolution_left[ist][op_i-1].empty())) continue;// this product matrix element is zero
    int ip = op_i;
    
    if (ip==0) nstate = praStates[ist];// we need to start from scratch because op_i iz zero
    else nstate = state_evolution_left[ist][ip-1];// we start from the previous operator and update from there on
    if (nstate.empty()) continue;
      
    bool empty=false;
    
    for (; ip<Operators.size(); ip++){        // go over the rest of the operators
      mstate.Evolve(nstate, cluster, Operators.exp_(ip)); // state evolution due to local Hamiltonian (is diagonal)
      int iop = Operators.typ(ip);
      nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);        // the next operator is applied
      if (nstate.empty()) {empty=true; break;}// maybe matrix element is zero
      state_evolution_left[ist][ip].SetEqual(nstate);      // it is not, store it for next time
    }

    if (!empty){ // from the last operator to beta
      mstate.Evolve(nstate, cluster, Operators.exp_last()); // The last evolution of the state due to H_loc.
      state_evolution_left[ist][ip++].SetEqual(mstate); // store it
    }
    if (ip>state_evolution_left.fullsize_Nd()) {cerr<<"Trying to resize to "<<ip<<" op_size="<<Operators.size()<<endl;}
    state_evolution_left[ist].resize(ip); // we have stored state evolution, remember the new size
    Number mm=0;
    if (!empty) mm = mstate.Project_to_Pra(praStates[ist], Projection[ist]);
    Trace[ist] = mm;
    ms += mm;
  }
  
  Number ms_right=0;
  for (int ist=0; ist<praStates.size(); ist++){// op_i is the earliest operator inserted - from this point on, the evolution needs to be changed
    int nsize = Operators.size();
    int ip = nsize-1-op_f;
    if (ip>0 && (state_evolution_right[ist].size()<ip || state_evolution_right[ist][ip-1].empty())) continue;// this product matrix element is zero
    if (ip==0) nstate = praStates[ist];// we need to start from scratch because op_i iz zero
    else nstate = state_evolution_right[ist][ip-1];// we start from the previous operator and update from there on
    if (nstate.empty()) continue;
    
    bool empty=false;
    
    if (!empty && ip<nsize){ // from the last operator to beta
      if (ip>0) mstate.Evolve(nstate, cluster, Operators.exp_(nsize-ip)); // The last evlution of the state due to H_loc.
      else 	mstate.Evolve(nstate, cluster, Operators.exp_last()); // The last evlution of the state due to H_loc.
      int iop = Operators.typ_transpose(nsize-1-ip);// HERE iop==-1 !!!!
      nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);        // the next operator is applied
      if (nstate.empty()) {empty=true;}// maybe matrix element is zero
      if (!empty) state_evolution_right[ist][ip++].SetEqual(nstate); // store it
    }
    if (!empty){
      for (; ip<Operators.size(); ip++){        // go over the rest of the operators
	mstate.Evolve(nstate, cluster, Operators.exp_(nsize-ip)); // state evolution due to local Hamiltonian (is diagonal)
	int iop = Operators.typ_transpose(nsize-1-ip);
	nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);        // the next operator is applied
	if (nstate.empty()) {empty=true; break;}// maybe matrix element is zero
	state_evolution_right[ist][ip].SetEqual(nstate);      // it is not, store it for next time
      }
    }
    if (!empty){
      mstate.Evolve(nstate, cluster, Operators.exp_(0)); // The last evlution of the state due to H_loc.
      state_evolution_right[ist][ip++].SetEqual(mstate); // store it
    }
    if (ip>state_evolution_right.fullsize_Nd()) {cerr<<"Trying to resize to "<<ip<<" op_size="<<Operators.size()<<endl;}
    state_evolution_right[ist].resize(ip); // we have stored state evolution, remember the new size
    Number mm=0;
    if (!empty) {
      mm = mstate.TraceProject(praStates[ist]); // need the trace. The last state needs to be projected to the first state
    }
    ms_right += mm;
  }

  if (fabs(fabs(divide(ms_right,ms))-1)>1e-6){
    cerr<<"ms left and ms right are not the same "<<ms<<" "<<ms_right<<endl;
    if (fabs(divide(ms_right,ms))>1) ms = ms_right;
  }
  
  if (common::PreciseP && cmp_P){ // computes probability more precisely
    double* __restrict__ Prob_ = Prob.MemPt();
    int npra = praStates.size();
    int plda = Prob.lda();

    int nsize = Operators.size();
    if (nsize==0){
      for (int ist=0; ist<npra; ist++)
	for (int m=0; m<Prob[ist].size(); m++)
	  Prob(ist,m) = cluster.P_atom(ist,m);
    }else{
      for (int ist=0; ist<npra; ist++)
	for (int m=0; m<Prob[ist].size(); m++)
	  Prob_[ist*plda+m] = 0.0;
      
      for (int ist=0; ist<npra; ist++){
	
	if (state_evolution_left[ist].size()<=Operators.size() || state_evolution_right[ist].size()<=Operators.size()) continue;
	if (op_i>0 && (state_evolution_left[ist].size()<op_i || state_evolution_left[ist][op_i-1].empty())) continue;// this product matrix element is zero

	for (int l=0; l<nsize-1; l++){
	  NState& lstate = state_evolution_left(ist,l);
	  NState& rstate = state_evolution_right(ist,nsize-l-2);
	  int istate = lstate.istate;
	  
	  if (lstate.istate != rstate.istate){
	    cout<<"ERROR: lstate.istate="<<lstate.istate<<" rstate.istate="<<rstate.istate<<endl;

	    cout<<"prastate="<<praStates[ist].istate<<endl;
	    for (int i=0; i<state_evolution_left[ist].size()-1; i++) {
	      cout<<"left "<<i<<" "<<state_evolution_left[ist][i].istate<<"    operator:"<<Operators.typ(i)<<endl;
	    }
	    for (int i=0; i<state_evolution_right[ist].size(); i++) {
	      cout<<"right "<<i<<" "<<state_evolution_right[ist][i].istate<<endl;
	    }
	  }
	  
	  double* __restrict__ lstateM = lstate.M.MemPt();
	  double* __restrict__ rstateM = rstate.M.MemPt();
	  int nk = lstate.M.size_Nd();
	  int ldal = lstate.M.lda();
	  int ldar = rstate.M.lda();
	  
	  //const function<double>& expn = Operators.exp_(l+1)[lstate.istate];
	  //function<double>& Prb = Prob[lstate.istate-1];
	  const double* __restrict__ expn = Operators.exp_(l+1)[istate].MemPt();
	  double* __restrict__ Prb = &Prob_[(istate-1)*plda];
	  double dtau = (Operators.t(l+1)-Operators.t(l))/common::beta;
	  double dtau_ms_mantisa = dtau/ms.mantisa;
	  double out_exponents = lstate.exponent+rstate.exponent-ms.exponent;
	  
	  for (int n=0; n<lstate.M.size_N(); n++){
	    double mantisa=0;
	    for (int m=0; m<nk; m++) mantisa += lstateM[n*ldal+m]*rstateM[n*ldar+m];
	    Prb[n] += mantisa*dtau_ms_mantisa*exp(expn[n]+out_exponents);
	  }
	}
    
	NState& lstate = state_evolution_left(ist,nsize);
	//function<double>& Prb = Prob[lstate.istate-1];
	double* __restrict__ Prb = &Prob_[(lstate.istate-1)*plda];
	double* __restrict__ lstateM = lstate.M.MemPt();
	int ldal = lstate.M.lda();
    
	double dtau = 1 - (Operators.t(nsize-1)-Operators.t(0))/common::beta;
	double dtau_ms_mantisa_exp = dtau/ms.mantisa * exp(lstate.exponent-ms.exponent);
	for (int n=0; n<lstate.M.size_N(); n++)
	  //Prb[n] += lstate.M(n,n)/ms.mantisa*exp(lstate.exponent-ms.exponent)*dtau;
	  Prb[n] += lstateM[n*ldal+n]*dtau_ms_mantisa_exp;
      }
    }
    double sum=0;
    for (int ist=0; ist<npra; ist++)
      for (int l=0; l<Prob[ist].size(); l++)
	sum += Prob_[ist*plda+l];
       //sum += Prob[ist][l];
  }else{
    double sum=0;
    for (int ist=0; ist<praStates.size(); ist++){
      for (int l=0; l<Prob[ist].size(); l++){
	Prob[ist][l] = divide(Projection[ist][l],ms);
	sum += Prob[ist][l];
      }
    }
  }
  
  return ms;
}

Number ComputeTrace(const function2D<NState>& state_evolution_left,
		    const function2D<NState>& state_evolution_right,
		    int op_i, int op_f,
		    const NOperators& Operators, const function1D<NState>& praStates,
		    const ClusterData& cluster,
		    int dsize)
{
  
  static NState nstate(common::max_size,common::max_size), mstate(common::max_size,common::max_size);
  Number ms_new=0;
  for (int ist=0; ist<praStates.size(); ist++){
    if (op_i>0 && (state_evolution_left[ist].size()<op_i || state_evolution_left[ist][op_i-1].empty())) continue;
    int ip = op_i;
    if (ip==0) nstate = praStates[ist];
    else nstate = state_evolution_left[ist][ip-1];
    if (nstate.empty()) continue;
      
    int nsize = Operators.size();

    // Go through states and check that matrix element is indeed nonzero
    // Do not multiply matrices yet.
    bool alive = true;
    int instate = nstate.istate;
    int ip0=ip;
    for (; ip0<=op_f; ip0++){
      int iop = Operators.tr_typ(ip0);
      instate = cluster.Fi(iop)[instate];
      if (instate==0) {	alive=false; break; }
    }
    if (!alive) continue;
    if (op_f<nsize+dsize-1){
      alive = (state_evolution_right[ist].size()>nsize+dsize-op_f-2 && instate==state_evolution_right[ist][nsize+dsize-op_f-2].istate);
    }else{
      alive = (praStates[ist].istate == instate);
    }
    if (!alive) continue;
        
    Number mm_new=0;
    bool empty=false;
    
    for (; ip<=op_f; ip++){
      mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
      int iop = Operators.tr_typ(ip);
      nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);
      if (nstate.empty()) {empty=true; break;}
    }
    if (!empty){
      if (op_f<nsize+dsize-1){
	if (state_evolution_right[ist].size()>nsize+dsize-op_f-2 &&
	    nstate.istate==state_evolution_right[ist][nsize+dsize-op_f-2].istate){
	  mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
	  mm_new = mstate.ScalarProduct(state_evolution_right[ist][nsize+dsize-op_f-2]);
	}else{
	  mm_new = 0;
	}
      }else{
	mstate.Evolve(nstate, cluster, Operators.tr_exp_last());
	mm_new = mstate.TraceProject(praStates[ist]);
      }
    }
    ms_new += mm_new;
    //cout<<"o_"<<ist<<" "<<mm_new<<endl;
  }
  return abs(ms_new);
}

Number ComputeGlobalTrace(const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, const function<int>& gflip_fl)
{
  static NState nstate, mstate;
  Number ms_new=0;
  int nsize = Operators.size();
  for (int ist=0; ist<praStates.size(); ist++){
    int ip = 0;
    nstate = praStates[ist];
    if (nstate.empty()) continue;
    
    Number mm_new=0;
    bool empty=false;
    
    for (; ip<nsize; ip++){
      mstate.Evolve(nstate, cluster, Operators.exp_(ip));
      int iop = Operators.typ(ip);
      int nfl = gflip_fl[iop/2];
      int nop = 2*nfl + iop%2;
      
      nstate.apply(cluster.FM(nop), cluster.Fi(nop), mstate);
      if (nstate.empty()) {empty=true; break;}
    }
    if (!empty){
      mstate.Evolve(nstate, cluster, Operators.exp_last());
      mm_new = mstate.TraceProject(praStates[ist]);
    }
    ms_new += mm_new;
  }
  return ms_new;
}

void ComputeMoments(const ClusterData& cluster, const mesh1D& omi, int iom_start, double nf,
		    const function1D<double>& mom, const function1D<double>& kaver, const function2D<double>& nt,
		    function2D<dcomplex>& Sigma, bool CorrectMoments=true, int Ncorrect=-1, int aom=3, double sderiv=0.1)
{
  if (Ncorrect<0) Ncorrect=cluster.N_unique_fl;
  stringstream script; script.precision(15);
  script<<"#!/usr/bin/perl\n";
  script<<"use Math::Trig;\nuse Math::Complex;\n\n"; 
  script<<"$nf="<<nf<<";\n";
  script<<"$T="<<(1/common::beta)<<";\n";
  script<<"$U="<<common::U<<";\n";
  script<<"$mu="<<common::mu<<";\n";
  script<<"$Nm="<<cluster.Nm<<";\n";
  script<<"$NK="<<cluster.N_unique_fl<<";\n";
  script<<"$NKc="<<Ncorrect<<";\n";
  for (int i=0; i<nt.size_N(); i++){
    script<<"@nt"<<i<<"=(";
    for (int j=0; j<nt.size_Nd()-1; j++) script<<nt[i][j]<<",";
    script<<nt[i].last()<<");\n";
  }
  for (int op=0; op<cluster.HF_M.size(); op++){
    script<<"@Op"<<op+1<<"=(";
    for (int ik=0; ik<cluster.N_unique_fl; ik++){
      script<<mom[ik+op*cluster.N_unique_fl];
      if (ik<cluster.N_unique_fl-1) script<<",";
      else script<<");\n";
    }
  }
  script<<"@epsk=(";
  for (int ik=0; ik<cluster.N_unique_fl; ik++){
    script<<cluster.epsk[ik];
    if (ik<cluster.N_unique_fl-1) script<<",";
    else script<<");\n";
  }
  
  script<<"@kaver=(";
  for (int ik=0; ik<kaver.size(); ik++){
    script<<kaver[ik];
    if (ik<kaver.size()-1) script<<",";
    else script<<");\n";
  }
  
  script<<"@om=(";
  for (int im=iom_start; im<omi.size(); im++){
    script<<omi[im];
    if (im<omi.size()-1) script<<",";
    else script<<");\n";
  }
  script<<"\n";
  script<<"for ($im=0; $im<=$#om; $im++){\n";
  script<<"  for ($k=0; $k<$NK; $k++){\n";
  script<<"    "<<cluster.RealSigma<<";"<<endl;
  script<<"    "<<cluster.ImagSigma<<";"<<endl;
  script<<"    $$Sig[$k][$im] = $ReSigma + $ImSigma*i;\n";
  script<<"  }\n";
  //  script<<"  print \"\\n\";\n";
  script<<"}\n";

  if (CorrectMoments){
    // Corrects high frequency by interpolation
    script<<endl<<endl;
    script<<"@sb=(";
    for (int ik=0; ik<cluster.N_unique_fl; ik++){
      dcomplex sb=0;
      for (int im=iom_start-aom; im<iom_start; im++) sb += Sigma[ik][im];
      sb/=aom;
      script<<sb.real()<<" ";
      if (sb.imag()>=0) script<<"+";
      script<<sb.imag()<<"*i";
      if (ik<cluster.N_unique_fl-1) script<<",";
      else script<<");"<<endl;
    }
    script<<endl;
    script<<"$omb = "<<omi[iom_start-aom/2]<<";"<<endl;
    script<<"$sderiv = "<<sderiv<<";"<<endl<<endl;
    script<<"for ($k=0; $k<$NKc; $k++) {$se[$k] = Re($$Sig[$k][$#om-1]);}"<<endl;
  
    script<<"for ($im=0; $im<=$#om; $im++){"<<endl;
    script<<"    for ($k=0; $k<$NKc; $k++){"<<endl;
    script<<"	$sr = $se[$k] + (Re($sb[$k])-$se[$k])*($omb/$om[$im])**2;"<<endl;
    script<<"	$$Sig[$k][$im] = $sr + Im($$Sig[$k][$im])*i;"<<endl;
    script<<"    }"<<endl;
    script<<"}"<<endl;
    script<<endl;
    script<<"for ($k=0; $k<$NKc; $k++){"<<endl;
    script<<"    for ($iw=1; $iw<$#om; $iw++){"<<endl;
    script<<"	$df0 = Im($$Sig[$k][$iw]-$sb[$k])/($om[$iw]-$omb);"<<endl;
    script<<"	$df1 = Im($$Sig[$k][$iw]-$$Sig[$k][$iw-1])/($om[$iw]-$om[$iw-1]);"<<endl;
    script<<"	if (abs($df0-$df1)<$sderiv) {last;}"<<endl;
    script<<"    }"<<endl;
    script<<"    $niw[$k] = $iw;"<<endl;
    script<<"    $se[$k] = Im($$Sig[$k][$iw]);"<<endl;
    script<<"    $oe[$k] = $om[$iw];"<<endl;
    script<<"}"<<endl;
    script<<endl;
    script<<"for ($k=0; $k<$NKc; $k++){"<<endl;
    script<<"    for ($im=0; $im<$niw[$k]; $im++){"<<endl;
    script<<"	$si = Im($sb[$k]) + ($se[$k]-Im($sb[$k]))*($om[$im]-$omb)/($oe[$k]-$omb);"<<endl;
    script<<"	$$Sig[$k][$im] = Re($$Sig[$k][$im]) + $si*i;"<<endl;
    script<<"    }"<<endl;
    script<<"}"<<endl;
  }
  script<<endl;
  script<<"for ($im=0; $im<=$#om; $im++){"<<endl;
  script<<"  print $om[$im], \"  \";"<<endl;
  script<<"  for ($k=0; $k<$NK; $k++){"<<endl;
  script<<"    print Re($$Sig[$k][$im]), \" \", Im($$Sig[$k][$im]), \"  \";"<<endl;
  script<<"  }\n";
  script<<"  print \"\\n\";\n";
  script<<"}\n";
  script<<endl;
  
  
  ofstream scrpt("pscript.pl");
  scrpt<<script.str()<<endl;
  
  string command;
  command = "perl pscript.pl > moments.temp";
  int ifail = system(command.c_str());

  ifstream inp("moments.temp");
  if (ifail){cerr<<"Seems your moments are not working or could not find perl!"<<endl;}
  for (int im=iom_start; im<omi.size(); im++){
    double omega;
    inp>>omega;
    if (fabs(omega-omi[im])>1e-6) {cerr<<"Did not succeed in reading moments. Exiting!"<<endl; break;}
    for (int ik=0; ik<cluster.N_unique_fl; ik++) inp>>Sigma(ik,im);
    inp.ignore(1000,'\n');
    if (!inp) {cerr<<"Did not succeed in reading moments. Exiting!"<<endl; break;}
  }
}




double ComputeExponentSum(int ist, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize)
{
  NState nstate = (op_i == 0) ? praStates[ist] : state_evolution_left[ist][op_i-1];
  int istate = nstate.istate;  // istate index superstates from 1 (vacuum = 0)
  double exp_sum = nstate.exponent;

  double mants=0;
  for (int i=0; i<nstate.M.size_N(); i++)
    for (int j=0; j<nstate.M.size_Nd(); j++)
      mants += fabs(nstate.M(i,j));
  exp_sum += log(mants);
  
  int last_op = Operators.size() + dsize - 1;  // index of last operator
  
  int ip = op_i;
  for ( ; ip<=op_f; ++ip) {
    const funProxy<double>& exps = Operators.tr_exp_(ip)[istate];
    exp_sum += *max_element(exps.begin(), exps.end());
    int iop = Operators.tr_typ(ip);
    istate = cluster.Fi(iop, istate);
  }

  if (op_f < last_op){
    const funProxy<double>& exps = Operators.tr_exp_(ip)[istate];
    exp_sum += *max_element(exps.begin(), exps.end());
    nstate = state_evolution_right[ist][last_op-op_f-1];
    exp_sum += nstate.exponent;

    double mants=0;
    for (int i=0; i<nstate.M.size_N(); i++)
      for (int j=0; j<nstate.M.size_Nd(); j++)
	mants += fabs(nstate.M(i,j));
    exp_sum += log(mants);
    
  }else{
    const funProxy<double>& exps = Operators.tr_exp_last()[istate];
    exp_sum += *max_element(exps.begin(), exps.end());
  }
  return exp_sum;
}

bool keep_superstate(int ist, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize)
{
  // simple integer lookup of whether evolution annihilates the superstate
  const int VACUUM = 0; // when superstate are indexed from 1, vacuum state is 0
  
  // check if this superstate has survived to point in time (op_i) where we want to insert first operator
  if (op_i>0 && (state_evolution_left[ist].size()<op_i || state_evolution_left[ist][op_i-1].empty())) return false;
  // superstate has survived to first operator insertion -- get evolution of praState up to this point
  NState nstate = (op_i == 0) ? praStates[ist] : state_evolution_left[ist][op_i-1];

  // check if state will survive until second operator insertion = end of modified time segment
  int evolved_state = nstate.istate;  // evolved_state index superstates from 1 (vacuum = 0)
  for (int ip = op_i; ip<=op_f; ++ip) {
    evolved_state = cluster.Fi(Operators.tr_typ(ip), evolved_state);
    if (evolved_state == VACUUM) return false;
  }

  int last_op = Operators.size() + dsize - 1;  // index of last operator
  if (op_f == last_op) {
    // if op_f is the last operator
    return (evolved_state == praStates[ist].istate);
  } else {
    // if op_f is not the last operator
    return (state_evolution_right[ist].size() > last_op-op_f-1 && evolved_state == state_evolution_right[ist][last_op-op_f-1].istate);// make sure we're not out of bounds
  }
}

bool compare_exp_sums(const pair<int, double> exp_sum1, const pair<int, double> exp_sum2)
{
  return (exp_sum1.second > exp_sum2.second);  // use > (instead of <) to place largest exponent at beginning of list
}

void ComputeExpTrace(list<pair<int,double> >& exp_sums, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize)
{
  // Algorithm
  // for istate = 1 .. N_superstates
  //   if superstate doesn't survive, continue
  //   otherwise, compute exponent sum and store
  // sort list of stored exponent sums
  // compute approximation to trace trace_approx = max(stored exponent sums)

  for (int ist=0; ist<praStates.size(); ist++){
    if ( keep_superstate(ist, state_evolution_left, state_evolution_right, op_i, op_f, Operators, praStates, cluster, dsize) ) {
      double exp_sum = ComputeExponentSum(ist, state_evolution_left, state_evolution_right, op_i, op_f, Operators, praStates, cluster, dsize);
      exp_sums.push_back( pair<int, double>(ist, exp_sum) );
      //LOGN0(TRACE, ist << "(" << exp_sum << ") ");
    }
  }
  //LOGN0(TRACE,endl);
  
  // sort list of stored exponent sums
  exp_sums.sort(compare_exp_sums); // largest element is first
}

Number ComputePartialTrace(int ist, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize)
{
  // compute trace for one string of superstates
  static NState nstate, mstate;
  Number mm_new = 0.0;

  // superstate has survived.  Get evolution of praState up to the point where we want to insert first operator
  int ip = op_i;
  nstate = (ip == 0) ? praStates[ist] : state_evolution_left[ist][ip-1];
  if (nstate.empty()) return mm_new;
  
  int nsize = Operators.size();
  int last_op = Operators.size() + dsize - 1;

  // evolve state with new operator
  for (; ip<=op_f; ip++){
    // first apply exponential (time) evolution
    mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
    // then apply creation/annihilation operator ("kink")
    int iop = Operators.tr_typ(ip);
    nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);
    if (nstate.empty()) return mm_new;
  }

  // apply final exponential (time) evolution
  if (op_f < last_op){
    if (state_evolution_right[ist].size()>last_op-op_f-1 && nstate.istate==state_evolution_right[ist][last_op-op_f-1].istate){
      mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
      mm_new = mstate.ScalarProduct(state_evolution_right[ist][last_op-op_f-1]);
    }
  }else{
    mstate.Evolve(nstate, cluster, Operators.tr_exp_last());
    mm_new = mstate.TraceProject(praStates[ist]);
  }
  return mm_new;
}


bool LazyComputeTrace(double P_r, const Number& matrix_element, double& ms, const list<pair<int,double> >& exp_sums, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int op_i, int op_f, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize)
{
  // Algorithm: see PRB , 005100 (2014).
  
  Number amatrix_element = abs(matrix_element);  // matrix element of the previous QMC step -- need for ratio
  
  Number approx_sum=0.0; // Sum of all approximate traces -- first we just add up exponents
  for (list<pair<int,double> >::const_iterator exp_sums_iter = exp_sums.begin(); exp_sums_iter!=exp_sums.end(); ++exp_sums_iter)
    approx_sum += Number(1,exp_sums_iter->second);
  
  double ms_max = divide(approx_sum,amatrix_element);  // this gives the first upper bound for probability
  if (ms_max<=P_r) {
    //if (P_exact>P_r) 
    //  cout<<"ERROR at the end a0 "<<P_exact<<" "<<P_r<<endl;
    ms=ms_max;
    return false; // Reject
  }
  
  
  Number total_trace=0.0;
  list<pair<int,double> >::const_iterator exp_sums_iter = exp_sums.begin();
  double ms_min;
  
  //int ii=0;
  while (exp_sums_iter!=exp_sums.end()){ 
    // Finds first non-zero partial trace and evaluate it exactly
    Number largest_partial_trace=0.0;
    while (exp_sums_iter != exp_sums.end() && largest_partial_trace.mantisa==0.0) {
      int ist = exp_sums_iter->first;
      Number partial_trace = ComputePartialTrace(ist, state_evolution_left, state_evolution_right, op_i, op_f, Operators, praStates, cluster, dsize);
      largest_partial_trace = partial_trace;
      ++exp_sums_iter;
    }
    total_trace += largest_partial_trace; // exact part of the trace up to now
    
    approx_sum=0.0; // Sum of all approximate traces -- the rest, which has not been computed exactly
    list<pair<int,double> >::const_iterator exp_sums_itern = exp_sums_iter;
    for (; exp_sums_itern!=exp_sums.end(); ++exp_sums_itern)
      approx_sum += Number(1,exp_sums_itern->second);
    
    ms_max = divide( abs(total_trace) + approx_sum, amatrix_element);  // upper bound for probability
    ms_min = divide( abs(total_trace) - approx_sum, amatrix_element);  // lower bound for probability
    //if (ms_max+1e-13<P_exact){ cout<<"ERROR 2b,"<<ii<<": "<<ms_max<<" "<<P_exact<<" "<<ms_max-P_exact<<endl;}
    //if (ms_min-1e-13>P_exact){ cout<<"ERROR 2c,"<<ii<<": "<<ms_min<<" "<<P_exact<<" "<<P_exact-ms_min<<endl;}
    //ii++;
    
    if (ms_max<=P_r) {
      //if (P_exact>P_r) 
      //  cout<<"ERROR at the end a "<<P_exact<<" "<<P_r<<endl;
      ms=ms_max;
      return false; // Reject
    }
    if (ms_min>=P_r) {
      //if (P_exact<P_r)
      //  cout<<"ERROR at the end b "<<P_exact<<" "<<P_r<<endl;
      ms=ms_min;
      return true;  //Accept
    }
  }
  if (fabs(ms_min-ms_max)>1e-13) cout<<"ERROR : Should be equal "<<ms_min<<" "<<ms_max<<endl;
  return P_r<ms_min;
  //return P_r<P_exact;
}

Number ComputeTraceFromScratch(const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster)
{
  static NState nstate, mstate;
  Number ms_new=0;
  int nsize = Operators.size();
  for (int ist=0; ist<praStates.size(); ist++){
    int ip = 0;
    nstate = praStates[ist];
    if (nstate.empty()) continue;
    
    Number mm_new=0;
    bool empty=false;
    
    for (; ip<nsize; ip++){
      mstate.Evolve(nstate, cluster, Operators.exp_(ip));
      int iop = Operators.typ(ip);
      nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);
      if (nstate.empty()) {empty=true; break;}
    }
    if (!empty){
      mstate.Evolve(nstate, cluster, Operators.exp_last());
      mm_new = mstate.TraceProject(praStates[ist]);
      //mm_new = mstate.Project_to_Pra(praStates[ist], Projection[ist]);
    }
    ms_new += mm_new;
  }
  return ms_new;
}








/*


bool compare_exp_sums_debug(const pair<int, double> exp_sum1, const pair<int, double> exp_sum2)
{
  return (exp_sum1.first < exp_sum2.first);  // use > (instead of <) to place largest exponent at beginning of list
}
Number ComputeTraceDEBUG(const list<pair<int,double> >& exp_sums, 
			 const function2D<NState>& state_evolution_left,
			 const function2D<NState>& state_evolution_right,
			 int op_i, int op_f,
			 const NOperators& Operators, const function1D<NState>& praStates,
			 const ClusterData& cluster,
			 int dsize)
{
  list<pair<int,double> > exp_sums_local(exp_sums);
  exp_sums_local.sort(compare_exp_sums_debug); // largest element is first
  
  for (list<pair<int,double> >::const_iterator ii=exp_sums_local.begin(); ii!=exp_sums_local.end(); ++ii){
    cout<<" "<< (ii->first) <<" "<< ii->second <<endl;
  }
  
  static NState nstate(common::max_size,common::max_size), mstate(common::max_size,common::max_size);
  Number ms_new=0;
  for (int ist=0; ist<praStates.size(); ist++){

    bool PRINT=false;
    if (ist==8) PRINT=true;


    if (PRINT){
      int ip = 0;
      nstate = praStates[ist];
      if (nstate.empty()) continue;
    
      Number mm_new=0;
      bool empty=false;
      for (; ip<Operators.size()+dsize; ip++){
	mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
	int iop = Operators.tr_typ(ip);
	nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);
	if (nstate.empty()) {empty=true; break;}
      }
      if (!empty){
	mstate.Evolve(nstate, cluster, Operators.tr_exp_last());
	mm_new = mstate.TraceProject(praStates[ist]);
      }
      cout<<"From scratch: "<<mm_new<<" "<<mm_new.balance()<<endl;
    }
    
    if (op_i>0 && (state_evolution_left[ist].size()<op_i || state_evolution_left[ist][op_i-1].empty())) continue;
    int ip = op_i;
    if (ip==0) nstate = praStates[ist];
    else nstate = state_evolution_left[ist][ip-1];
    if (nstate.empty()) continue;
      
    int nsize = Operators.size();

    // Go through states and check that matrix element is indeed nonzero
    // Do not multiply matrices yet.
    bool alive = true;
    int instate = nstate.istate;
    int ip0=ip;
    for (; ip0<=op_f; ip0++){
      int iop = Operators.tr_typ(ip0);
      instate = cluster.Fi(iop)[instate];
      if (instate==0) {	alive=false; break; }
    }
    if (!alive) continue;
    if (op_f<nsize+dsize-1){
      alive = (state_evolution_right[ist].size()>nsize+dsize-op_f-2 && instate==state_evolution_right[ist][nsize+dsize-op_f-2].istate);
    }else{
      alive = (praStates[ist].istate == instate);
    }
    if (!alive) continue;
    
    Number mm_new=0;
    bool empty=false;

    if (PRINT){
      cout<<"Start:"<<nstate.istate<<" "<<nstate.exponent<<endl;
    }
    for (; ip<=op_f; ip++){
      mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
      int iop = Operators.tr_typ(ip);
      nstate.apply(cluster.FM(iop), cluster.Fi(iop), mstate);
      if (nstate.empty()) {empty=true; break;}
    }
    if (!empty){
      if (op_f<nsize+dsize-1){
	if (state_evolution_right[ist].size()>nsize+dsize-op_f-2 &&
	    nstate.istate==state_evolution_right[ist][nsize+dsize-op_f-2].istate){
	  mstate.Evolve(nstate, cluster, Operators.tr_exp_(ip));
	  mm_new = mstate.ScalarProduct(state_evolution_right[ist][nsize+dsize-op_f-2]);
	}else{
	  mm_new = 0;
	}
      }else{
	mstate.Evolve(nstate, cluster, Operators.tr_exp_last());
	mm_new = mstate.TraceProject(praStates[ist]);
      }
    }
    ms_new += mm_new;
    cout<<"ist="<<ist<<" "<<mm_new<<" "<<mm_new.balance()<<endl;
  }

  exit(0);
  return abs(ms_new);
}
*/
