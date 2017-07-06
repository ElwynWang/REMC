void SetupMatrixStuff(void);
void Fold(void);

void SetupMatrixStuff(void) {
  int i, j, k;
  struct cell *the_cell;

  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      for(k=0; k<MATRIX_SIZE; k++) 
	the_matrix[i][j][k].natoms = 0;

  for(i=0; i<natoms; i++) { 
    FindLatticeCoordinates(&native[i]);

    the_cell = &the_matrix[native[i].X][native[i].Y][native[i].Z];
    the_cell->atom_list[(the_cell->natoms)++] = i;

  }
}

void Fold(void) {
  int i;

  int irep;
  int itmp;
  float etmp;

  char rmsd_filename[100];
  char temp_filename[100];
  float E_pot_now=0.0;
  float E_tor_now=0.0;
  float E_sct_now=0.0;
  float E_aro_now=0.0;
  float E_hbond_now=0.0;

  /*rmsd by jyang*/
  struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE];
  struct alignment align;

  NFRAG = 1;
  align.seqptr[1].x1=1;
  align.structptr[1].x1=1;
  align.seqptr[1].x2=nresidues;
  align.structptr[1].x2=nresidues;

  mc_flags.init=!NO_NEW_CLASHES;
  for(i=0; i<natoms; i++)
    CopyAtom(native[i],&orig_native[i]);

  CenterProtein(&native,natoms);  
  SetupMatrixStuff();  

  Contacts();  
  fprintf(STATUS,"INITIAL CLASHES\t%d\n",nclashes); 
  if (!mc_flags.init) {
    fprintf(STATUS,"\nturning off clashes:\n");
    TurnOffNativeClashes(1);  // init.h
  }
  fprintf(STATUS,"\n");

  /* setup initial energies */
  //align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &native_rms);
  for(i=0;i<nalign;++i)
   {
    struct_f1[i+1].CA.x=struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;   // the crystal structure of the protein
    struct_f1[i+1].CA.y=struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
    struct_f1[i+1].CA.z=struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
    struct_f2[i+1].CA.x=native[native_residue[map_to_seq[i]].CA].xyz.x;        // the model structure of the protein
    struct_f2[i+1].CA.y=native[native_residue[map_to_seq[i]].CA].xyz.y;
    struct_f2[i+1].CA.z=native[native_residue[map_to_seq[i]].CA].xyz.z;
   }
  native_rms = getrms(struct_f1,struct_f2,align);

  ResetEnergies(0);
  rms_RMSDmin = 100.;
  mcstep_RMSDmin = 0;
  mcstep_Emin = 0;
  Emin = prev_E;
  Emin_pot = prev_E_pot;
  Emin_hbond = prev_E_hbond;
  Emin_tor = prev_E_tor;
  Emin_sct = prev_E_sct;
  Emin_aro = prev_E_aro;
  fprintf(STATUS,"ENERGY = %.2f(clashes) + %.2f(rmsd) + %.2f(potential) + %.2f(Aro) + %.2f(hbond) + %.2f(tor) + %.2f(sct)\n\n",
    weight_clash, weight_rms, weight_potential, ARO_WEIGHT, weight_hbond, TOR_WEIGHT, SCT_WEIGHT);
  for(i=0; i<natoms; i++)
   {
    prev_native[i] = native[i];
    native_Emin[i] = prev_native[i];
   }
  for (i = 0; i < nresidues; ++i)
    cur_rotamers[i] = 0;

  naccepted=0;
  nrejected=0;
  nothers=0;
  n_sidechain_accepted=0;

  fprintf(STATUS,"         step #    energy contact rmsd     potnl    sctor      hbond       Aro   torsion accept reject others   temp   frag_acc yang_acc loop_acc local_acc frag_step yang_step loop_step local_step\n---------------------------------------------------------------------------------------------\n");

  /* main folding loop */

  for(mcstep = 0; mcstep < MC_STEPS; mcstep++) {
    // Replica Exchange MC 
    if (mcstep%MC_REPLICA_STEPS == 0) 
     {

        for(i=0;i<nprocs;i++) replica_index[i]=i;
        ierr=MPI_Allgather(&E,1,MPI_FLOAT,Enode,1,MPI_FLOAT,mpi_world_comm);

        if(myrank == 0) {
          for(irep=0;irep<MAX_EXCHANGE;irep++){
            sel_num= (int) (drand48()*(nprocs-2));  // drand48 generates a even distributed number between 0 and 1. Why nprocs-2 rather than nprocs-1

            //fprintf(STATUS,"irep : %d, sel_num : %d\n", irep, sel_num);
            //fflush(STATUS);

            delta_E = Enode[sel_num+1]-Enode[sel_num];  // delta energy 
            delta_T = 1.0/Tnode[sel_num+1]-1.0/Tnode[sel_num]; // delta temperature
            delta_all = delta_E * delta_T;
            if(delta_all >= 0 || drand48() < expf(delta_all)){
              itmp=replica_index[sel_num];
              replica_index[sel_num]=replica_index[sel_num+1];
              replica_index[sel_num+1]=itmp;
              etmp=Enode[sel_num];
              Enode[sel_num]=Enode[sel_num+1];
              Enode[sel_num+1]=etmp;
              accepted_replica[sel_num]++;
            } else {
              rejected_replica[sel_num]++;
            }
          }
        }

        ierr=MPI_Bcast(replica_index,nprocs,MPI_INT,0,mpi_world_comm);
        ierr=MPI_Bcast(Enode,nprocs,MPI_FLOAT,0,mpi_world_comm);
        ierr=MPI_Bcast(accepted_replica,nprocs,MPI_INT,0,mpi_world_comm);
        ierr=MPI_Bcast(rejected_replica,nprocs,MPI_INT,0,mpi_world_comm);

        //fprintf(STATUS,"Replica Index\n");
        //for(i=0;i<nprocs;i++) fprintf(STATUS,"%5d %5d %8.3f\n", i, replica_index[i], Enode[i]);

        fprintf(STATUS,"RPLC %10ld E : %8.3f FROM %2d(%5.3f) E : %8.3f, accepted : %5d, rejected : %5d\n", mcstep, E, replica_index[myrank], Tnode[replica_index[myrank]], Enode[myrank], accepted_replica[myrank], rejected_replica[myrank]);
         fflush(STATUS);

        for(i=0;i<natoms;i++) {
          buf_out[3*i]=native[i].xyz.x;
          buf_out[3*i+1]=native[i].xyz.y;
          buf_out[3*i+2]=native[i].xyz.z;
        }

        for(i=0;i<nprocs;i++){
          if(replica_index[i] != i) {
            if(myrank == i) {
              ierr=MPI_Recv(buf_in,3*natoms,MPI_FLOAT,replica_index[i],(i+2),mpi_world_comm,&mpi_status);
              //fprintf(STATUS,"Receiving coord. to %2d nodes..., ierr : %d\n", replica_index[i],ierr);
              //fprintf(STATUS,"MPIRcv: %8.3f %8.3f\n", buf_in[0], buf_in[natoms-1]);
            } else if (myrank == replica_index[i]) {
              ierr=MPI_Send(buf_out,3*natoms,MPI_FLOAT,i,(i+2),mpi_world_comm);
              //fprintf(STATUS,"Sending coord. to %2d nodes..., ierr : %d\n", i,ierr);
              //fprintf(STATUS,"MPISnd: %8.3f %8.3f\n", buf_out[0], buf_out[natoms-1]);
            }
          }
        }

        if(replica_index[myrank] != myrank) {
          for(i=0;i<natoms;i++) {
            native[i].xyz.x=buf_in[3*i];
            native[i].xyz.y=buf_in[3*i+1];
            native[i].xyz.z=buf_in[3*i+2];
          }
        }
        //fprintf(STATUS,"Finished replica exchange...\n");
      
      /* Re-center */
//      if (!USE_ROTAMERS) {
	CenterProtein(&native,natoms);
	SetupMatrixStuff();
	for(i=0; i<natoms; i++)
	  prev_native[i] = native[i];
//      }
      /* Re-set values that are susceptible to round-off error */
      Contacts();
      if (!mc_flags.init)
	TurnOffNativeClashes(0);
      ResetEnergies(0);
      GetChi();
     }

    if (backbone_accepted ==1)
     {
      for(i=0;i<nalign;++i)
       {
        struct_f1[i+1].CA.x=struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;
        struct_f1[i+1].CA.y=struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
        struct_f1[i+1].CA.z=struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
        struct_f2[i+1].CA.x=native[native_residue[map_to_seq[i]].CA].xyz.x;
        struct_f2[i+1].CA.y=native[native_residue[map_to_seq[i]].CA].xyz.y;
        struct_f2[i+1].CA.z=native[native_residue[map_to_seq[i]].CA].xyz.z;
       }
      native_rms = getrms(struct_f1,struct_f2,align);
     }
      
    /* Print Output */
    if (mcstep%MC_PRINT_STEPS==0) {
      //align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &native_rms);
    /*rmsd by jyang*/
      for(i=0;i<nalign;++i)
       {
        struct_f1[i+1].CA.x=struct_native[struct_residue[map_to_struct[i]].CA].xyz.x;
        struct_f1[i+1].CA.y=struct_native[struct_residue[map_to_struct[i]].CA].xyz.y;
        struct_f1[i+1].CA.z=struct_native[struct_residue[map_to_struct[i]].CA].xyz.z;
        struct_f2[i+1].CA.x=native[native_residue[map_to_seq[i]].CA].xyz.x;
        struct_f2[i+1].CA.y=native[native_residue[map_to_seq[i]].CA].xyz.y;
        struct_f2[i+1].CA.z=native[native_residue[map_to_seq[i]].CA].xyz.z;
       }
      native_rms = getrms(struct_f1,struct_f2,align);
      TypeContacts();
      E_pot_now = FullAtomEnergy();
      E_sct_now = sctenergy();
      E_tor_now = torsionenergy();
      E_aro_now = aromaticenergy();
      E_hbond_now = HydrogenBonds();
      fprintf(STATUS,"STEP %10ld  %8.2f %6d %5.2f %9.2f %9.2f %9.2f    %6.2f %8.2f %6.2f %6.2f %6.2f %6.3f   %6.2f %6.2f %6.2f %6.2f %6d %6d %6d %6d\n", 
	  mcstep, E, ncontacts, native_rms, E_pot_now, E_sct_now, E_hbond_now, E_aro_now, E_tor_now,
	  100*(Float)naccepted/(Float)MC_PRINT_STEPS, 100*(Float)nrejected/(Float)MC_PRINT_STEPS, 100*(Float)nothers/(Float)MC_PRINT_STEPS, MC_TEMP,\
    100*(Float)naccepted_frag/(Float)n_frag,100*(Float)naccepted_yang/(Float)n_yang,100*(Float)naccepted_loop/(Float)n_loop,100*(Float)naccepted_local/(Float)n_local,n_frag,n_yang,n_loop,n_local); 
      fflush(STATUS);
      n_sidechain_accepted=0;
      naccepted = 0;
      nrejected=0;
      nothers=0;
      naccepted_local=0;
      naccepted_loop=0;
      naccepted_yang=0;
      naccepted_frag=0;
      n_frag=0;
      n_yang=0;
      n_local=0;
      n_loop=0;
    } 

    if (mcstep%10000==0) {
      /* Re-center */
      CenterProtein(&native,natoms);
      SetupMatrixStuff();
      for(i=0; i<natoms; i++)
        prev_native[i] = native[i];

      /* Re-set values that are susceptible to round-off error */
      Contacts();
      if (!mc_flags.init)
	TurnOffNativeClashes(0);
      ResetEnergies(mcstep);
      GetChi();
    }
    
    /* Output Structure */
    if (PRINT_PDB) {
      if (mcstep%MC_PDB_PRINT_STEPS==0 && myrank==0) {
	sprintf(temp_filename,"%s.%ld",pdb_out_file,mcstep);  
    //    if ((MC_TEMP < 0.18) || (mcstep > MC_STEPS - 100000))
	  PrintPDB(temp_filename);
      }
    }

    if (E<Emin)
     {
      mcstep_Emin = mcstep;
      Emin = E;
      Emin_pot = E_pot;
      Emin_hbond = E_hbond;
      Emin_tor = E_tor;
      Emin_aro = E_aro;
      rms_Emin = native_rms;
      for(i=0;i<natoms;i++)
	native_Emin[i] = native[i];
     }

    if (native_rms < rms_RMSDmin)
     {
      mcstep_RMSDmin = mcstep;
      rms_RMSDmin = native_rms;
      E_RMSDmin = E;
      E_RMSDmin_pot = E_pot;
      E_RMSDmin_hbond = E_hbond;
      E_RMSDmin_tor = E_tor;
      E_RMSDmin_sct = E_sct;
      E_RMSDmin_aro = E_aro;
      for(i=0;i<natoms;i++)
	native_RMSDmin[i] = native[i];
     }
    /* Make a move */
    MakeMove(STEP_SIZE,USE_GLOBAL_BB_MOVES);
  }
  if (myrank==0){
  sprintf(temp_filename,"%s_Emin.pdb",pdb_out_file);
  PrintPDB_Emin(temp_filename);}
  fprintf(STATUS,"\nMC step at Emin: %10ld\n", mcstep_Emin);
  fprintf(STATUS,"Emin:%8.2f  Emin_pot:%8.2f  E_hb:%7.2f  E_tor:%7.2f  E_sct:%7.2f E_aro:%7.2f ",
          Emin, Emin_pot,  Emin_hbond, Emin_tor, Emin_sct, Emin_aro);
  fprintf(STATUS,"rmsd at Emin:%8.2f  ", rms_Emin);
  fprintf(STATUS,"Pdb file at Emin: %s\n", temp_filename);
  if (myrank==0){
  sprintf(rmsd_filename,"%s_RMSDmin.pdb",pdb_out_file);
  PrintPDB_RMSDmin(rmsd_filename);}
  fprintf(STATUS,"\nMC step at RMSDmin: %10ld\n", mcstep_RMSDmin);
  fprintf(STATUS,"rms_Rmin:%8.2f  E_RMSDmin:%8.2f E_Rmin_pot:%8.2f E_Rhb:%8.2f E_Rtor:%8.2f E_Rsct:%8.2f E_Raro:%8.2f\n",
         rms_RMSDmin, E_RMSDmin, E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro);
  fprintf(STATUS,"Pdb file at RMSDmin: %s\n", rmsd_filename);

  return;  
}
