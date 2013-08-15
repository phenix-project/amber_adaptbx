#ifndef _thisheaderfile_H_
#define _thisheaderfile_H_


extern "C"
{

        int getmdgxfrc(const char *tpname, const char *crdname, const double *PhenixCoords, \
                           double* target, double * gradients);
}


#endif
