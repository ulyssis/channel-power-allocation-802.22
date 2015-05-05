#include<mex.h>
#include <stdlib.h>
#include <stdio.h>
//int system(const char *string);
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
int callScript()
{
system("cd /home/li/work/tools/lindo/lindoapi/samples/c/dica && ./whitecat_compoundcoefficient" );
return 0;
}


void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[]){
callScript();
}
