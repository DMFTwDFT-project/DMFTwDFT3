
void ClusterData::Read_for_HB1(ifstream& gout, double mu, double U)
{
  cnsize=0;
  if (!QHB1) return;

  string str;
  clog<<"Continuing with HB1 after "<<Osize<<" operators read!"<<endl;
  int wNm, wN_ifl, wmax_size;
  gout>>wNm>>cnsize>>wN_ifl>>wmax_size;
  if (wNm!=Nm) {cerr<<"Nm="<<wNm<<" for HB1 and original Nm="<<Nm<<" are not the same"<<endl;exit(1);}
  if (wN_ifl!=N_ifl) {cerr<<"Nm="<<wN_ifl<<" for HB1 and original Nm="<<N_ifl<<" are not the same"<<endl;exit(1);}
  getline(gout,str); 
  getline(gout,str); // # ind   N   K   Jz size
    
    
  cNs.resize(cnsize+1);
  cmsize_.resize(cnsize+1);
  cindx.resize(cnsize+1);
  cF_i.resize(N_flavors,cnsize+1);
  cEne.resize(cnsize+1);
  cF_M.resize(N_flavors,cnsize+1);
  cF_i=0;
  
  for (int i=1; i<=cnsize; i++){
    int it1, in, ik, isize, tindx;
    double dsz;
    gout>>it1>>tindx>>in>>ik>>dsz>>isize;
    if (it1!=i){cerr<<"Something wrong in parsing cix file. it1!=i!"<<endl;}
    cNs[i] = in;
    cmsize_[i] = isize;
    cindx[i] = tindx;
    for (int ib=0; ib<N_flavors; ib++){
      int ist;
      gout>>ist;
      if (ist!=0){
	cF_i(ib,i)=ist;
      }
    }
    cEne[i].resize(isize);
    for (int is=0; is<isize; is++) gout>>cEne[i][is];
    double spin;
    for (int is=0; is<isize; is++) gout>>spin;
    gout.ignore(1000,'\n');
  }
  getline(gout,str); // # matrix elements
    
  for (int i=1; i<=cnsize; i++){
    for (int ib=0; ib<N_flavors; ib++){
      int it, jt, size1, size2;
      gout>>it>>jt>>size1>>size2;
      if (it!=i || jt!=cF_i(ib,i)) cerr<<"Something wrong reading cix file"<<endl;
      cF_M(ib,i).resize(size2,size1); // constructor changes i->ii
      //      int ii = cF_i(ib,i);
      for (int i1=0; i1<size1; i1++){
	for (int i2=0; i2<size2; i2++){
	  double m;
	  gout>>m;
	  cF_M(ib,i)(i2,i1)=m; // constructor only
	}
      }
      gout.ignore(1000,'\n');
    }
  }
  for (int i=1; i<=cnsize; i++){
    double dE = -cNs[i]*mu + 0.5*cNs[i]*(cNs[i]-1)*U;
    for (size_t m=0; m<cEne[i].size(); m++) cEne[i][m] += dE;
  }

}

void ClusterData::HB1(double beta, double mu, const mesh1D& iom, const function2D<double>& AProb,
		      int nom, int aom, const function2D<dcomplex>& Delta, function2D<dcomplex>& Sigma, const function1D<bool>& nonzero, const function1D<bool>& nonzerou, double sderiv=0.1, int nom_small=300)
{
  // Interpolate probabilities to the large atomic base
  vector<function1D<double> > Prob(cnsize);
  for (int j=0; j<cnsize; j++){
    int isize = cmsize_[j+1];
    Prob[j].resize(isize);
    Prob[j]=0;
    int indx = cindx[j+1];
    if (indx==0) continue;
    for (int m=0; m<msize(indx); m++){
      Prob[j][m] = AProb[indx-1][m];
    }
  }

  // Create a small logarithimic mesh of 30 points
  mesh1D ioms;
  CreateLogMesh(1, nom_small, 1./beta, iom, ioms);

  // Computes atomic Green's function
  function2D<dcomplex> Gh(N_unique_fl,ioms.size());
  Gh=0;
  function1D<dcomplex> gh(ioms.size()); 
  for (int ifl=0; ifl<N_ifl; ifl++){
    if (! nonzero[ifl]) continue; // Some GF are projected out. Should be ignored
    //    clog<<ifl<<" ";
    for (int i1=0; i1<ifl_dim[ifl]; i1++){
      for (int i2=0; i2<ifl_dim[ifl]; i2++){
	int ia = fl_from_ifl(ifl,i1);
	int ib = fl_from_ifl(ifl,i2);
	int tfl = tfl_index[ifl][i1][i2];
	int bfl = bfl_index[ifl][tfl];
	int sgn = sign[ifl][tfl];
	int dcmp = conjg[ifl][tfl];
	gh=0;
	for (int i=1; i<=cnsize; i++){
	  int j = cF_i(ia,i);
	  int jb = cF_i(ib,i);
	  if (j!=jb || j==0) continue;
	  for (int im=0; im<cmsize_[i]; im++){
	    double Pm = Prob[i-1][im];
	    for (int jm=0; jm<cmsize_[j]; jm++){
	      double Pn = Prob[j-1][jm];
	      double mm = cF_M(ia,i)(jm,im)*cF_M(ib,i)(jm,im)*(Pn+Pm);
	      double dE = cEne[j][jm]-cEne[i][im];
	      dcomplex ci = dcomplex(0,1);
	      double dE2 = sqr(dE);
	      double mdE = mm*dE;
	      for (int io=0; io<ioms.size(); io++){
		double ome = ioms[io];
		double den = 1/(ome*ome+dE2);
		gh[io].real() -= mdE*den;// gh = mm/(iom-dE)
		gh[io].imag() -= mm*ome*den;
	      }
	    }
	  }
	}
	for (int io=0; io<ioms.size(); io++){
	  dcomplex f = sgn*gh[io];
	  if (dcmp) f = f.conj();
	  Gh[bfl][io] += f;
	}
      }
    }
  }
  for (int ifl=0; ifl<N_unique_fl; ifl++) Gh[ifl] *= (1./fl_deg[ifl]);
  //  clog<<endl;

  // Computes Hubbard I self-energy
  function2D<dcomplex> Gf(N_unique_fl,ioms.size());
  function2D<dcomplex> Gf_1(N_unique_fl,ioms.size());
  Gf_1=0;
  vector<spline1D<dcomplex> > Sigh(N_unique_fl);
  for (size_t fl=0; fl<Sigh.size(); fl++){
    Sigh[fl].resize(ioms.size());
    for (int im=0; im<ioms.size(); im++) Sigh[fl][im]=0;
  }
  
  for (int im=0; im<ioms.size(); im++){
    dcomplex iomega(0,ioms[im]);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) = Gh(fl, im);
    for (int ifl=0; ifl<N_ifl; ifl++){
      if (nonzero[ifl]) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf, Gf_1);
    }
    for (int fl=0; fl<N_unique_fl; fl++) Gf_1(fl,im) /= fl_deg[fl];
    for (int fl=0; fl<N_unique_fl; fl++){
      if (nonzerou[fl]) Sigh[fl][im] = (iomega+mu)*Id[fl]-epsk[fl]-Gf_1(fl,im);
    }
  }
  
  // brisi!
  Gf=0;
  Gf_1=0;
  ofstream hout("g_hb1.dat");
  ofstream sout("s_hb1.dat");
  ofstream rout("g_hb0.dat");
  for (int im=0; im<ioms.size(); im++){
    dcomplex iomega(0,ioms[im]);
    for (int fl=0; fl<N_unique_fl; fl++) Gf_1(fl,im) = (iomega+mu)*Id[fl]-epsk[fl]-Sigh[fl][im]-Delta[fl](iom.Interp(ioms[im]));
    for (int ifl=0; ifl<N_ifl; ifl++) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf_1, Gf);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) /= fl_deg[fl];
    hout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) hout<<Gf[fl][im]<<" ";
    hout<<endl;
    sout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) sout<<Sigh[fl][im]<<" ";
    sout<<endl;
    rout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) rout<<Gh[fl][im]<<" ";
    rout<<endl;
  }
  //brisi!

  
  for (size_t fl=0; fl<Sigh.size(); fl++){
    int n = ioms.size()-1;
    dcomplex df0 = (Sigh[fl][1]-Sigh[fl][0])/(ioms[1]-ioms[0]);
    dcomplex dfn = (Sigh[fl][n]-Sigh[fl][n-1])/(ioms[n]-ioms[n-1]);
    Sigh[fl].splineIt(ioms, df0, dfn);
  }

  // Interpolates between low energy QMC and high energy HBI self-energy
  for (int ifl=0; ifl<N_unique_fl; ifl++){
    if (!nonzerou[ifl]) continue;
    dcomplex sb=0;
    for (int im=nom-aom; im<nom; im++) sb += Sigma[ifl][im];
    sb/=aom;
    double ob = iom[nom-aom/2];
    // Imaginary part
    int is=nom;
    for (; is<iom.size()-1; is++){
      intpar p0 = ioms.Interp(iom[is]);
      intpar p1 = ioms.Interp(iom[is-1]);
      double df0 = (Sigh[ifl](p0)-sb).imag()/(iom[is]-ob);
      double df1 = (Sigh[ifl](p0)-Sigh[ifl](p1)).imag()/(iom[is]-iom[is-1]);
      if (fabs(df0-df1)<sderiv) break;
    }
    double se = Sigh[ifl](ioms.Interp(iom[is])).imag();
    double oe = iom[is];
    for (int im=nom; im<is; im++) Sigma[ifl][im].imag() = sb.imag() + (se-sb.imag())*(iom[im]-ob)/(oe-ob);
    for (int im=is; im<iom.size(); im++) Sigma[ifl][im].imag() = Sigh[ifl](ioms.Interp(iom[im])).imag();
    // Real part
    double sinf = Sigh[ifl].last().real();
    for (int im=nom; im<iom.size(); im++) Sigma[ifl][im].real() = sinf + sqr(ob/iom[im])*(sb.real()-sinf);
  }

  
  Gf.resize(N_unique_fl,iom.size());
  Gf_1.resize(N_unique_fl,iom.size());
  Gf=0;
  Gf_1=0;
  for (int im=0; im<iom.size(); im++){
    dcomplex iomega(0,iom[im]);
    for (int fl=0; fl<N_unique_fl; fl++)
      if (nonzerou[fl]) Gf_1(fl,im) = (iomega+mu)*Id[fl]-epsk[fl]-Sigma(fl,im)-Delta(fl,im);
    for (int ifl=0; ifl<N_ifl; ifl++)
      if (nonzero[ifl]) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf_1, Gf);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) /= fl_deg[fl];
  }
  
  ofstream tout("g_qmc.dat");
  for (int io=0; io<iom.size(); io++){
    tout<<iom[io]<<" ";
    for (int ifl=0; ifl<N_unique_fl; ifl++) tout<<Gf[ifl][io]<<" ";
    tout<<endl;
  }
}
