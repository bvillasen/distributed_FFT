#ifndef GLOBAL_H 
#define GLOBAL_H 


#if PRECISION == 1
#ifndef FLOAT_TYPEDEF_DEFINED
typedef float Real;
#endif //FLOAT_TYPEDEF_DEFINED
#endif //PRECISION == 1
#if PRECISION == 2
#ifndef FLOAT_TYPEDEF_DEFINED
typedef double Real;
#endif //FLOAT_TYPEDEF_DEFINED
#endif //PRECISION == 2




#endif //GLOBAL_H
