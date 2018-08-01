#ifndef __SUPPORT_HH
#define __SUPPORT_HH

#include <stdlib.h>

#include <iostream>
#include <string>

template<typename T1, typename T2>
void assert_equal ( const T1 &x , const T2 &y , std::string context ) 
{
  if ( x != y ) {
    std::cout << "Error while " << context << std::endl;
    std::cout << x << " != " << y << std::endl;
    exit ( 1 );
  }
}

#endif //#ifndef __SUPPORT_HH
