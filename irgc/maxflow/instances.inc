#include "graph_irgc.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <captype, tcaptype, flowtype>
// IMPORTANT: 
//    flowtype should be 'larger' than tcaptype 
//    tcaptype should be 'larger' than captype

template class GraphIrgc<short,short,int>;
template class GraphIrgc<short,short,long>;
template class GraphIrgc<float,float,double>;

template class GraphIrgc<int,int,int>;
template class GraphIrgc<short,int,int>;
template class GraphIrgc<float,float,float>;
template class GraphIrgc<double,double,double>;

