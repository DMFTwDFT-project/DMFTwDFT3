#ifdef _MPI
#include <mpi.h>
using namespace std;

void Reduce(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<dcomplex>& Gd, function2D<dcomplex>& Sd,
	    function2D<double>& AverageProbability, double& asign, function1D<double>& nlc, function1D<double>& kaver, function2D<double>& susc,
	    function2D<double>& Gtau, function5D<dcomplex>& VertexH, function5D<dcomplex>& VertexF, function1D<int>& Gd_deg,
	    bool cmp_vertex, bool QHB2)
{
  function2D<double> cAverageProbability;
  function1D<double> chistogram;
  function1D<double> cnlc;
  function1D<double> ckaver;
  function2D<dcomplex> cGd;
  function2D<dcomplex> cSd;
  function2D<double> cSusc;
  function2D<double> cGtau;
  function1D<int> cGd_deg;
  double casign;
  if (my_rank==Master){
    chistogram.resize(histogram.size());
    cnlc.resize(nlc.size());
    cAverageProbability.resize(AverageProbability.fullsize_N(),AverageProbability.fullsize_Nd());
    cGd.resize(Gd.fullsize_N(),Gd.fullsize_Nd());
    if (QHB2) cSd.resize(Sd.fullsize_N(),Sd.fullsize_Nd());
    cGtau.resize(Gtau.fullsize_N(),Gtau.fullsize_Nd());
    ckaver.resize(kaver.size());
    cSusc.resize(susc.fullsize_N(), susc.fullsize_Nd());
    cGd_deg.resize(Gd_deg.size());
  }
  
  MPI::COMM_WORLD.Reduce(histogram.MemPt(), chistogram.MemPt(), histogram.size(), MPI_DOUBLE, MPI_SUM, Master);

  MPI::COMM_WORLD.Reduce(AverageProbability.MemPt(), cAverageProbability.MemPt(), AverageProbability.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);

  MPI::COMM_WORLD.Reduce(&asign, &casign, 1, MPI_DOUBLE, MPI_SUM, Master);
  
  MPI::COMM_WORLD.Reduce(nlc.MemPt(), cnlc.MemPt(), nlc.size(), MPI_DOUBLE, MPI_SUM, Master);
  
  MPI::COMM_WORLD.Reduce(kaver.MemPt(), ckaver.MemPt(), kaver.size(), MPI_DOUBLE, MPI_SUM, Master);
  
  MPI::COMM_WORLD.Reduce(Gtau.MemPt(), cGtau.MemPt(), Gtau.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);

  MPI::COMM_WORLD.Reduce(Gd.MemPt(), cGd.MemPt(), Gd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master);

  if (QHB2) MPI::COMM_WORLD.Reduce(Sd.MemPt(), cSd.MemPt(), Sd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master);
  
  MPI::COMM_WORLD.Reduce(susc.MemPt(), cSusc.MemPt(), susc.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);

  MPI::COMM_WORLD.Reduce(Gd_deg.MemPt(), cGd_deg.MemPt(), Gd_deg.size(), MPI_INT, MPI_SUM, Master);
  
  if (cmp_vertex){
    function2D<dcomplex> cVertex(VertexH.N3, VertexH.N4);
    int psize = VertexH.N3*VertexH.N4;
    
    for (int i0=0; i0<VertexH.N0; i0++){
      for (int i1=0; i1<VertexH.N1; i1++){
	for (int i2=0; i2<VertexH.N2; i2++){
	  
	  dcomplex* f = &VertexH(i0,i1,i2,0,0);
	  MPI::COMM_WORLD.Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master);
	  
	  for (int i3=0; i3<VertexH.N1; i3++)
	    for (int i4=0; i4<VertexH.N2; i4++)
	      VertexH(i0,i1,i2,i3,i4) = cVertex(i3,i4)*(1./mpi_size);
	  
	  f = &VertexF(i0,i1,i2,0,0);
	  MPI::COMM_WORLD.Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master);
	  
	  for (int i3=0; i3<VertexH.N1; i3++)
	    for (int i4=0; i4<VertexH.N2; i4++)
	      VertexF(i0,i1,i2,i3,i4) = cVertex(i3,i4)*(1./mpi_size);

	  
	}
      }
    }
  }
  
  if (my_rank==Master){
    histogram = chistogram;
    histogram *= (1./mpi_size);
    AverageProbability = cAverageProbability;
    asign = casign;
    AverageProbability *= (1./mpi_size);
    nlc = cnlc;
    nlc *= (1./mpi_size);
    kaver = ckaver;
    kaver *= (1./mpi_size);
    Gd = cGd;
    if (QHB2) Sd = cSd;
    Gtau = cGtau;
    Gtau *= (1./mpi_size);
    susc = cSusc;
    susc *= (1./mpi_size);
    asign *= (1./mpi_size);
    Gd_deg = cGd_deg;
  }
}

void MPI_Init(int argc, char* argv[], int& my_rank, int& mpi_size, int& Master)
{
  MPI::Init(argc, argv);
  my_rank = MPI::COMM_WORLD.Get_rank();
  mpi_size = MPI::COMM_WORLD.Get_size();
  Master = 0;
  
  std::cout << "Hello World! I am " << my_rank << " of " << mpi_size << std::endl;
}

void MPI_finalize()
{MPI::Finalize();}

#else

using namespace std;

void Reduce(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<dcomplex>& Gd, function2D<dcomplex>& Sd,
	    function2D<double>& AverageProbability, double& asign, function1D<double>& nlc,
	    function1D<double>& kaver, function2D<double>& susc, function2D<double>& Gtau,
	    function5D<dcomplex>& VertexH, function5D<dcomplex>& VertexF, function1D<int>& Gd_deg, bool cmp_vertex, bool QHB2){}

void MPI_Init(int argc, char* argv[], int& my_rank, int& mpi_size, int& Master)
{
  my_rank = 0;
  mpi_size = 1;
  Master = 0;
  std::cout<<"Not parallel!"<<std::endl;
}

void MPI_finalize(){}
#endif
