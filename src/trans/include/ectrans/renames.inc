!#define PASTER(x,y) x ## _ ## y
!#define EVALUATOR(x,y)  PASTER(x,y)
#define PASTE2(x,y) x/**/y
#define EVALUATOR(x,y) PASTE2(PASTE2(x,_),y)
#define SUFFIXIZE(fun) EVALUATOR(fun,SYMBOLSUFFIX)
#define upwind3d_gpubc SUFFIXIZE(upwind3d_gpubc)
