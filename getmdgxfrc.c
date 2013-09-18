#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdgx.h"
#include "time.h"


trajcon CreateTrajCon(){
	trajcon tj; //trajectory control data (input file params)
	InitBasicTrajcon(&tj);
	return tj;
}

uform LoadTopology(char *tpname, trajcon* tj){
    uform U;     /*topology, direct and recip space controls and lookup tables
                   convolution support */
	U = InitPotential(tpname, 8.0, tj);
	return U;
} 

mdsys CreateMDSys(char *crdname, uform* U){
	mdsys thisMD; /*MD structs: coords, cell grid, energy decomp, timings*/
	thisMD.crd = (crdname[0] != '\0') ?
		ReadRst(&U->tp, crdname) : ReadRst(&U->tp, "inpcrd");
	return thisMD;	
}

 
void LoadPhenixCoordToGrid(char* crdname, uform *U, trajcon *tj, 
                       const double *PhenixCoords, mdsys* thisMDptr )
{
  /*** Replace old coords with new Phenix coords ***/
  int k;
  for (k=0; k<(*thisMDptr).crd.natom; k++) {
        (*thisMDptr).crd.loc[k*3]=PhenixCoords[k*3];
        (*thisMDptr).crd.loc[k*3+1]=PhenixCoords[k*3+1];
        (*thisMDptr).crd.loc[k*3+2]=PhenixCoords[k*3+2];
  }
  ImageBondedGroups(&(*thisMDptr).crd, &U->tp);
  InitHistory(&(*thisMDptr).crd);
  /*** Create the cell grid and prepare reciprocal space support ***/
  (*thisMDptr).CG = CreateCellGrid(&(*thisMDptr).crd, &U->dcinp, &U->rcinp, &U->tp, tj, 0);
  PrepPME(&(*thisMDptr).CG, &U->rcinp, &(*thisMDptr).crd);
  U->PPk = CreateBCKit(&U->rcinp, &U->rcinp.QL[0], &(*thisMDptr).crd, &U->tp,
                       FFTW_ESTIMATE);
#ifdef MPI
  LinkCellGrid(&(*thisMDptr).CG, &(*thisMDptr).crd, &U->rcinp);
#else
  LinkCellGrid(&(*thisMDptr).CG, &U->rcinp);
#endif
  /*** Load the cell grid ***/
  AtomsToCells(&(*thisMDptr).crd, &(*thisMDptr).CG, &U->tp);
}


int getmdgxfrc(char *tpname, char *crdname, const double *PhenixCoords,
               double* target, double * gradients, uform* Uptr, 
               trajcon* tjptr, mdsys* MDptr)
{

  LoadPhenixCoordToGrid(crdname, Uptr, tjptr, PhenixCoords, MDptr);
  //Compute forces
  InitExecon( &(*MDptr).etimers );
  MMForceEnergy(Uptr, MDptr, tjptr);

  //Populate gradients array	
  int i, j;
  cell *C;
  for (i = 0; i < (*MDptr).CG.ncell; i++) {
    C = &(*MDptr).CG.data[i];
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
  
  //Populate energy components vector
  *target=(*MDptr).sysUV.etot;
  *(target+1)=(*MDptr).sysUV.bond;
  *(target+2)=(*MDptr).sysUV.angl;
  *(target+3)=(*MDptr).sysUV.dihe;
  *(target+4)=(*MDptr).sysUV.elec;
  *(target+5)=(*MDptr).sysUV.vdw12+(*MDptr).sysUV.vdw6;
  *(target+6)=(*Uptr).tp.withH.nbond;
  *(target+7)=(*Uptr).tp.withH.nangl;
  *(target+8)=(*Uptr).tp.withH.ndihe;
  
  return 0;
}
