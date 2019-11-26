#include <iostream>
#include <cmath>

class Number{
public:
  double mantisa;
  double exponent;
  Number(double mantisa_=0, double exponent_=0) : mantisa(mantisa_), exponent(exponent_){};
  Number(const Number& z){ mantisa=z.mantisa; exponent=z.exponent;}
  Number& operator=(const Number& z){ mantisa=z.mantisa; exponent=z.exponent; return *this;}
  double dbl() const {return mantisa*exp(exponent);}
  double exp_dbl() const { return log(abs(mantisa))+exponent;}
  Number& operator+= (const Number& z);
  Number& operator-= (const Number& z);
  Number& operator*= (const Number& z);
  Number& operator/= (const Number& z);
  
  Number& balance();
  
  friend Number operator* (const Number& x, const Number& y) {return Number(x.mantisa*y.mantisa, x.exponent+y.exponent);};
  friend Number operator* (const Number& x, double y) {return Number(x.mantisa*y, x.exponent);};
  friend Number operator+ (const Number& x, const Number& y);
  friend Number operator- (const Number& x, const Number& y);
  friend Number operator/ (const Number& x, const Number& y) {return Number(x.mantisa/y.mantisa, x.exponent-y.exponent);};
  friend double divide(const Number& x, const Number& y) {return x.mantisa/y.mantisa*exp(x.exponent-y.exponent);};
  
  friend bool operator == (const Number& x, const Number& y) {return x.mantisa==y.mantisa && x.exponent==y.exponent;};
  friend bool operator == (const Number& x, double y) {return x.mantisa*exp(x.exponent) == y;}
  friend bool operator == (double x, const Number& y) {return y==x;};
  
  friend bool operator != (double x, const Number& y) {return !(x==y);}
  friend bool operator != (const Number& x, double y) {return !(x==y);}
  friend bool operator != (const Number& x, const Number& y) { return !(x==y);}
  friend std::ostream& operator<< (std::ostream& stream, const Number& z) { stream<<z.mantisa<<" "<<z.exponent<<" "; return stream;};
  friend Number sqrt(const Number& z) { return Number(sqrt(z.mantisa), z.exponent/2);};
  friend bool isnan(const Number& z) { return ::isnan(z.mantisa) || ::isnan(z.exponent);}
  friend Number abs(const Number& z) { return Number(fabs(z.mantisa),z.exponent);}
  friend double log(const Number& z) { return log(z.mantisa)+z.exponent;};
  friend Number  operator- (const Number& z) { return Number(-z.mantisa,z.exponent);};

  friend bool operator > (const Number& x, const Number& y) {
    return (x.exponent > y.exponent) ? (x.mantisa > y.mantisa * exp(y.exponent-x.exponent)) : (x.mantisa * exp(x.exponent-y.exponent) > y.mantisa);
  }
  friend bool operator < (const Number& x, const Number& y) {
    return (x.exponent > y.exponent) ? (x.mantisa < y.mantisa * exp(y.exponent-x.exponent)) : (x.mantisa * exp(x.exponent-y.exponent) < y.mantisa);
  }
  friend bool operator >= (const Number& x, const Number& y) { return !operator<(x, y); }
  friend bool operator <= (const Number& x, const Number& y) { return !operator>(x, y); }
  
};

Number& Number::operator+= (const Number& z)
{
  if (z.exponent<exponent)
    mantisa += z.mantisa*exp(z.exponent-exponent);
  else {
    mantisa = z.mantisa+mantisa*exp(exponent-z.exponent);
    exponent=z.exponent;
  }
  return *this;
}
Number& Number::operator-= (const Number& z)
{
  if (z.exponent<exponent)
    mantisa -= z.mantisa*exp(z.exponent-exponent);
  else{
    mantisa = -z.mantisa+mantisa*exp(exponent-z.exponent);
    exponent=z.exponent;
  }
  return *this;
} 
Number& Number::operator*= (const Number& z)
{
  mantisa *= z.mantisa;
  exponent += z.exponent;
  return *this;
}
Number& Number::operator/= (const Number& z)
{
  mantisa /= z.mantisa;
  exponent -= z.exponent;
  return *this;
}
Number& Number::balance()
{
  if (mantisa==0) {exponent=0; return *this;}
  exponent += log(fabs(mantisa));
  mantisa = (mantisa>0) ? 1 : -1;
  return *this;
}

Number operator+ (const Number& x, const Number& y)
{ return (x.exponent>y.exponent) ? Number(x.mantisa+y.mantisa*exp(y.exponent-x.exponent), x.exponent) : Number(y.mantisa+x.mantisa*exp(x.exponent-y.exponent), y.exponent);}
Number operator- (const Number& x, const Number& y)
{ return (x.exponent>y.exponent) ? Number(x.mantisa-y.mantisa*exp(y.exponent-x.exponent), x.exponent) : Number(-y.mantisa+x.mantisa*exp(x.exponent-y.exponent), y.exponent);}

Number balance(const Number& x)
{
  if (x.mantisa == 0) return Number(0.0, 0.0);
  return Number( (x.mantisa > 0) ? 1. : -1., x.exponent + log(fabs(x.mantisa)) );
}

