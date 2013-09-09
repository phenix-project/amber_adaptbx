#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdgx.h"
 
mdsys LoadCoordToGrid1(char* crdname, uform *U, trajcon *tj, const double *PhenixCoords)
{
  mdsys thisMD;
  /*** Read coordinates from disk ***/
  thisMD.crd = (crdname[0] != '\0') ?
    ReadRst(&U->tp, crdname) : ReadRst(&U->tp, "inpcrd");
  int k;
  for (k=0; k<thisMD.crd.natom; k++) {
        thisMD.crd.loc[k*3]=PhenixCoords[k*3];
        thisMD.crd.loc[k*3+1]=PhenixCoords[k*3+1];
        thisMD.crd.loc[k*3+2]=PhenixCoords[k*3+2];
  }
  ImageBondedGroups(&thisMD.crd, &U->tp);
  InitHistory(&thisMD.crd);
  /*** Create the cell grid and prepare reciprocal space support ***/
  thisMD.CG = CreateCellGrid(&thisMD.crd, &U->dcinp, &U->rcinp, &U->tp, tj, 0);
  PrepPME(&thisMD.CG, &U->rcinp, &thisMD.crd);
  U->PPk = CreateBCKit(&U->rcinp, &U->rcinp.QL[0], &thisMD.crd, &U->tp,
                       FFTW_ESTIMATE);
#ifdef MPI
  LinkCellGrid(&thisMD.CG, &thisMD.crd, &U->rcinp);
#else
  LinkCellGrid(&thisMD.CG, &U->rcinp);
#endif
  /*** Load the cell grid ***/
  AtomsToCells(&thisMD.crd, &thisMD.CG, &U->tp);
  return thisMD;
}

int getmdgxfrc(char *tpname, char *crdname, const double *PhenixCoords,
                           double* target, double * gradients)
{
  trajcon tj;  /*trajectory control data (input file params)
                                        (topology file struct is prmtop)*/
  uform U;     /*potential energy function data (also uses prmtop)
                                 topology, direct and recip space controls and lookup tables
                                 convolution support */
  mdsys MD;    /*MD structs: coords, cell grid, energy decomp, timings*/

  InitBasicTrajcon(&tj);
  U = InitPotential(tpname, 8.0, &tj);
  MD = LoadCoordToGrid1(crdname, &U, &tj, PhenixCoords);
  // Compute forces
  InitExecon(&MD.etimers);

  //NOTE:a failure occurs here presumably when phenix passes sites_cart that
  // had different numer of atoms than the amber crd file...
  //printf("HELLO 1\n");
  MMForceEnergy(&U, &MD, &tj);
  //printf("HELLO 2\n");

  int i, j;
  cell *C;
  for (i = 0; i < MD.CG.ncell; i++) {
    C = &MD.CG.data[i];
    for (j = 0; j < C->nr[0]; j++) {
      //~ printf("Atom %4d at [ %9.4lf %9.4lf %9.4lf ] with force "
             //~ "[ %9.4lf %9.4lf %9.4lf ]\n", C->data[j].id, C->data[j].loc[0],
             //~ C->data[j].loc[1], C->data[j].loc[2], C->data[j].frc[0],
             //~ C->data[j].frc[1], C->data[j].frc[2]);
          gradients[C->data[j].id*3]=C->data[j].frc[0];
          gradients[C->data[j].id*3+1]=C->data[j].frc[1];
          gradients[C->data[j].id*3+2]=C->data[j].frc[2];
    }
  }

  //~ printf("Total energy = %12.6lf from:\n", MD.sysUV.etot);
  //~ printf("Bonds     = %12.6lf    Elec  = %12.6lf\n", MD.sysUV.bond,
         //~ MD.sysUV.elec);
  //~ printf("Angles    = %12.6lf    vdW   = %12.6lf\n", MD.sysUV.angl,
         //~ MD.sysUV.vdw12 + MD.sysUV.vdw6);
  //~ printf("Dihedral  = %12.6lf\n", MD.sysUV.dihe);
  //~ printf("%12.6f\n", MD.sysUV.eptot);
  *target=MD.sysUV.etot;
  *(target+1)=MD.sysUV.bond;
  *(target+2)=MD.sysUV.angl;
  *(target+3)=MD.sysUV.dihe;
  *(target+4)=MD.sysUV.elec;
  *(target+5)=MD.sysUV.vdw12+MD.sysUV.vdw6;
  return 0;
}
