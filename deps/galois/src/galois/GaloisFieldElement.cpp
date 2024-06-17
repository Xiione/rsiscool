#include "GaloisFieldElement.h"

namespace galois
{

   GaloisFieldElement::GaloisFieldElement(GFSymbol v)
   {
      poly_value = v & GF.size();
   }

   GaloisFieldElement::GaloisFieldElement(const GaloisFieldElement& gfe)
   {
      poly_value  = gfe.poly_value;
   }


   std::ostream& operator << (std::ostream& os, const GaloisFieldElement& gfe)
   {
      os << gfe.poly_value;
      return os;
   }


   GaloisFieldElement operator+(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result += b;
      return result;
   }


   GaloisFieldElement operator-(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result -= b;
      return result;
   }


   GaloisFieldElement operator*(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result *= b;
      return result;
   }


   GaloisFieldElement operator*(const GaloisFieldElement& a, const GFSymbol& b)
   {
      GaloisFieldElement result  = a;
      result *= b;
      return result;
   }


   GaloisFieldElement operator*(const GFSymbol& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = b;
      result *= a;
      return result;
   }


   GaloisFieldElement operator/(const GaloisFieldElement& a, const GaloisFieldElement& b)
   {
      GaloisFieldElement result  = a;
      result /= b;
      return result;
   }


   GaloisFieldElement operator^(const GaloisFieldElement& a, const int& b)
   {
      GaloisFieldElement result  = a;
      result ^= b;
      return result;
   }

}
