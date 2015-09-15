#include "IRGC.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <nodeid, labelid, captype>
// IMPORTANT: 
//    nodeid must be integer type, can be signed or unsigned
//	  labelid must be a signed integer type 
//    captype must be a signed data type

//template class IRGC<int,short,int>;
//template class IRGC<int,int,int>;
//template class IRGC<unsigned int,short,int>;
//template class IRGC<unsigned int,unsigned short,int>;
//template class IRGC<unsigned int,int,int>;
//template class IRGC<unsigned long,int,int>;
//template class IRGC<int,short,float>;
//template class IRGC<unsigned int,int,float>;
//template class IRGC<unsigned long,int,double>;

template class IRGC<long,int,float>;