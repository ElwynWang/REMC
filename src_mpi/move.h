void BackboneMove(Float);
void LoopBackboneMove(Float);
void LocalBackboneMove(Float);
void SidechainMove();
void MakeMove(Float, Float);

/*==================================================*/
/*                  main program                    */
/*==================================================*/
//double bef_phi[MAXSEQUENCE], bef_psi[MAXSEQUENCE];
void MakeMove(Float step_size, Float use_global_bb_moves) {
  
  int reject, del, N, M, i;
  sidechain_step = 0;
  int use_yang = 0;
  int use_frag =0;
  int use_loop=0;
  int use_local=0;
  int n_soln = 0;  // ?
//  float hb_before_move = 0.;
//  float hb_after_move = 0.;
//  int j;
//  float e;
/*rmsd by jyang*/
  struct backbone struct_f2[MAXSEQUENCE];
  struct alignment align;

  NFRAG = 1;
  align.seqptr[1].x1=1;
  align.structptr[1].x1=1;
  align.seqptr[1].x2=nresidues;
  align.structptr[1].x2=nresidues;
  backbone_accepted = 0;
 
  do 
  {

    reject=0;
    mc_flags.clashed=0;
    nomove=0;
    sidemovedone = 0;
    
    /* make backbone move */  // LoopBackboneMove(0.33333),YANG_move(0.5),others:LocalBackboneMove
    if (sidechain_step == 0){

      if ((use_global_bb_moves) && (use_global_bb_moves  > 0.001))
      {
	       fprintf(STATUS,"USE_GLOBAL_BB_MOVES is turuned on!!!");
	       exit(1);
      }

      fragmove_weight=exp(fragmove_accepted/(fragmove_accepted+fragmove_rejected));
      
      if (drand48()<FRAG_MOVE)
      {
        if (!USE_CRANK_MOVE)
          FragMove();
        else
          CrankFragMove();
        use_yang=0;
        use_frag=1;
        use_loop=0;
        use_local=0;
        n_frag+=1;
      }
      
      else if (drand48()<CLUSTER_MOVE)  // default=0.33333 in define.h
      {
        //fprintf(STATUS, "LoopBackbondMove():\n");
        LoopBackboneMove(step_size);
        use_yang = 0;
        use_frag=0;
        use_loop=1;
        use_local=0;
        n_loop+=1;
      }

      else
      {
        if (YANG_MOVE)  // default=0.5 in cfg
        {
  	      if (drand48()<YANG_MOVE)  
	        {
            //fprintf(STATUS, "integloop():\n");
            integloop(step_size, &n_soln);  // n_soln: number of alternative loop closure solutions.
	          use_yang = 1;
            use_frag=0;
            use_loop=0;
            use_local=0;
            n_yang+=1;
	        }

	        else
	        {
             //fprintf(STATUS, "LocalBackboneMove():\n");
             LocalBackboneMove(step_size);
	           use_yang = 0;
             use_frag=0;
             use_loop=0;
             use_local=1;
             n_local+=1;
	        }
        }

        else
        {
          //fprintf(STATUS, "LocalBackboneMove():\n");
          LocalBackboneMove(step_size);
          use_yang = 0;
          use_frag=0;
          use_loop=0;
          use_local=1;
          n_local+=1;
        }
      }
    }

    else 
    {
      if(total_ntorsions!=0)
      {
        //fprintf(STATUS, "SidechainMove():\n");
        SidechainMove();
      }
      use_yang = 0;
    }
    
    if (!mc_flags.init && mc_flags.clashed)
    {
        reject = 1;
        if ((nomove==0)&&(sidechain_step == 0))
          nrejected++;
    }
    else if ((use_yang == 1) && (n_soln==0))
    {
        reject = 1;
        if (sidechain_step == 0)
          nothers++;
    }
    else 
    { 
        delta_nclashes=0;
        if (mc_flags.init)  // three lines for mc_flags_init
          for(i=0; i<total_pairs2; i++) 
            delta_nclashes+=data[cd[i].a][cd[i].b].delta_clashes-data[cd[i].a][cd[i].b].clashes; //??

        dE=0; dE_pot = 0; dE_hbond = 0; dE_tor = 0; dE_aro = 0; dE_sct = 0;  // if this step can be accepted, then calculate the delta energy

        if (sidechain_step==0)
          dE_tor = torsionenergy() - prev_E_tor;

        if((sidechain_step!=0)&& sidemovedone ) 
        {
          dE_sct = sctenergy() - prev_E_sct;
        }

        dE_aro = aromaticenergy() - prev_E_aro;
        delta_contacts=0;
        for(i=0; i < total_pairs; i++) 
        {
          N = ab[i].a; M = ab[i].b;
	       del = data[N][M].delta_contacts-data[N][M].contacts;
	       delta_contacts += del;
	       dE_pot+=((Float) del)*potential[native[N].smogtype][native[M].smogtype];
        }
        
        if (weight_hbond) 
        {
          if (sidechain_step==0)
	           dE_hbond = FoldHydrogenBonds() - prev_E_hbond;
        }

        if ((weight_rms!=0.0)&&(sidechain_step==0))
        {
            for(i=0;i<nalign;++i)
            {
	             struct_f2[i+1].CA.x=native[native_residue[map_to_seq[i]].CA].xyz.x;
               struct_f2[i+1].CA.y=native[native_residue[map_to_seq[i]].CA].xyz.y;
               struct_f2[i+1].CA.z=native[native_residue[map_to_seq[i]].CA].xyz.z;
            }
        }
        //	align_drms(native, native_residue, struct_native, struct_residue, map_to_seq, map_to_struct, nalign, &bb_rms);
        dE = weight_potential*dE_pot + weight_clash*delta_nclashes + weight_hbond*dE_hbond + TOR_WEIGHT*dE_tor + ARO_WEIGHT*dE_aro + SCT_WEIGHT*dE_sct;
       
        /* Metropolis */
        if (dE <= 0 || drand48() < exp(-(dE)/MC_TEMP)) 
        {
          Update();
          if ((use_frag==1) && (sidechain_step)==0)
            fragmove_accepted+=1.0;
          
      
          //ResetEnergies();
          if ((nomove==0)&&(sidechain_step == 0))
	        {
	          naccepted++;
	          backbone_accepted = 1;
            if (use_frag==1)
              naccepted_frag+=1;
            else if (use_yang==1)
              naccepted_yang+=1;
            else if (use_loop==1)
              naccepted_loop+=1;
            else if (use_local==1)
              naccepted_local+=1;
	        }
	        else 
	          n_sidechain_accepted++;
            //    fprintf(STATUS,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", E, 0.2*E_pot+0.8*E_hbond, E_pot, E_hbond, dE, dE_pot, dE_hbond);
        }
        else
        {
	        reject=1;
          if ((nomove==0)&&(sidechain_step == 0))
	          nrejected++;
        }
    }
  
    if (reject) 
    { 
      Restore();
      if ((use_frag==1) && (sidechain_step)==0)
        fragmove_rejected+=1.0;
      //fprintf(STATUS, "rejected ...\n");
      //ResetEnergies();
    }

    for(i=0;i<all_rotated_natoms;i++)
    {
      is_rotated[all_rotated_atoms[i]]=0;
    }

    /*test
    check_phipsi();
    
    for (i=0;i<nresidues-2;i++){
      if (bef_phi[i]!=cur_phi[i])
        printf ("bef_phi:%8.3f cur_phi:%8.3f\n",bef_phi[i],cur_phi[i]);
      if (bef_psi[i]!=cur_psi[i])
        printf ("bef_psi:%8.3f cur_psi:%8.3f\n",bef_psi[i],cur_psi[i]);}
    
    for (i=0;i<nresidues-2;i++)
    {
      if (abs(cur_phi[i]*deg2rad-native_residue[i+1].phi)>0.001 || (native_residue[i+1].psi-cur_psi[i]*deg2rad)>0.001){
        printf ("Index:%d,cur_phi:%8.3f, native_phi:%8.3f; cur_psi:%8.3f, native_psi:%8.3f; cur_phi-native_phi:%8.3f, cur_psi-native_psi:%8.3f; Uneuqal Error!!!\n",\
          i,cur_phi[i]*deg2rad,native_residue[i+1].phi,cur_psi[i]*deg2rad,native_residue[i+1].psi,cur_phi[i]*deg2rad-native_residue[i+1].phi,cur_psi[i]*deg2rad-native_residue[i+1].psi);
      }
      else{
        printf ("Index:%d,cur_phi:%8.3f, native_phi:%8.3f; cur_psi:%8.3f, native_psi:%8.3f; cur_phi-native_phi:%8.3f, cur_psi-native_psi:%8.3f;\n",\
          i,cur_phi[i]*deg2rad,native_residue[i+1].phi,cur_psi[i]*deg2rad,native_residue[i+1].psi,cur_phi[i]*deg2rad-native_residue[i+1].phi,cur_psi[i]*deg2rad-native_residue[i+1].psi);
      }
    }*/
  } 

  while(sidechain_step++ < SIDECHAIN_MOVES);  // one backbone move and one sidechain move in each step
  
  return;
  
}

/*==================================================*/


void LocalBackboneMove(Float step_size) {  // rotate phi or psi for a single residue, then psi or psi values of the shorter part of the whole chain will be changed according to the rotation axis of the selected residue; so it's a global move; should check line731-794 in init.h
  int i;

  mc.loop_size = 1;
  all_rotated_natoms = 0;
  total_pairs=total_pairs2=0;
  total_hbond_pairs=0;
  nomove = 0;

  mc.sel_res_num = (int)(drand48()*nresidues);
  if(is_template[mc.sel_res_num]==1)
   {
    nomove = 1;
    nothers++;
    return;
   }
  if (native_residue[mc.sel_res_num].amino_num==14) { /* cannot rotate around proline phi */
    if (mc.sel_res_num != nresidues-1) 
      mc.is_phi = 0; 
    else /* if proline is last residue, no rotation possible */
     {
      nomove = 1;
      nothers++;
      return;
     }
  }
  else if (mc.sel_res_num == 0)
    mc.is_phi = 0;
  else if (mc.sel_res_num == nresidues-1)
    mc.is_phi = 1;
  else    mc.is_phi = (int) (drand48()*2);
  //fprintf(STATUS,"mainchain move at %5d\n", mc.sel_res_num);

  step_size = step_size*GaussianNum();;  // GaussianNum in misc_util.h : to generate a number but don't understand the value,maybe gaussian distribution

  BackboneMove(step_size);
  
  UpdateLattice(rotate_natoms[mc.is_phi][mc.sel_res_num], rotate_atom[mc.is_phi][mc.sel_res_num]);

  for(i=0;i<rotate_natoms[mc.is_phi][mc.sel_res_num];i++){
    is_rotated[rotate_atom[mc.is_phi][mc.sel_res_num][i]]=1;
  }
  
  NewDeltaContacts(rotate_natoms[mc.is_phi][mc.sel_res_num], rotate_atom[mc.is_phi][mc.sel_res_num], not_rotated[mc.is_phi][mc.sel_res_num]);
  return;

}

void LoopBackboneMove(Float absolute_step_size) { // it's the fragment-based move set; like the LocalBackboneMove, it's also a global move; the only difference from LocalBackBoneMove is to sleect the changing angle range by cluster rather than random range from Gaussaion distribution; should check line796-944 in init.h

  int a,b,c,d, i, j;
  Float step_phi, step_psi;
  Float desire_phi, desire_psi;
  int cluster_bin;
  float use_cluster;

  /* every 1 moves, a new triplet of residues is selected. */
  /* each move, the phi/psi angles of a single residue are changed */
 
  mc.loop_size = 1; 

  /* triplet is selected */
  if (0 == 0) {
    mc.sel_triplet = (int) (drand48()*TOTAL_TRIPLE_LOOP_MOVES)+TOTAL_SINGLE_LOOP_MOVES+TOTAL_DOUBLE_LOOP_MOVES;// select a triplet number
 
    if (residue_triplets[mc.sel_triplet].a > nresidues/2.0) {
      mc.selected[0] = residue_triplets[mc.sel_triplet].a;
      mc.selected[1] = residue_triplets[mc.sel_triplet].b;
      mc.selected[2] = residue_triplets[mc.sel_triplet].c;
    }
    else {   
      if (residue_triplets[mc.sel_triplet].c>=0) {
	mc.selected[0] = residue_triplets[mc.sel_triplet].c;
	mc.selected[1] = residue_triplets[mc.sel_triplet].b;
	mc.selected[2] = residue_triplets[mc.sel_triplet].a;
      }
      else if (residue_triplets[mc.sel_triplet].b>=0) {
	mc.selected[0] = residue_triplets[mc.sel_triplet].b;
	mc.selected[1] = residue_triplets[mc.sel_triplet].a;
	mc.selected[2] = residue_triplets[mc.sel_triplet].c;
      } 
      else {
	mc.selected[0] = residue_triplets[mc.sel_triplet].a;
	mc.selected[1] = residue_triplets[mc.sel_triplet].b;
	mc.selected[2] = residue_triplets[mc.sel_triplet].c;
      }   
    }
  }
  if(is_template[mc.selected[0]]==1)
   {
    nomove = 1;
    nothers++;
    return;
   }      

  all_rotated_natoms = 0;
  total_pairs=total_pairs2=0;
  total_hbond_pairs=0;

  /* move-cycle defines the residue to be changed */
  i = move_cycle%1;  // important,in backbone.h =0
  move_cycle = i;
//  ++move_cycle;
  
  check_phipsi();
      
  cluster_bin = (int)(NOCLUSTERS*drand48()); // 30 in define.h
//  mc.selected[i] = 92;
//  cluster_bin = 1;
  desire_phi = cluster_phi[native_residue[mc.selected[i]].amino_num][cluster_bin];
  desire_psi = cluster_psi[native_residue[mc.selected[i]].amino_num][cluster_bin];
//fprintf(STATUS,"%d\n", cluster_bin);
  /* rotates phi and psi */
  if(secstr[mc.selected[i]]=='H')
    use_cluster = 0.5;
  else if(secstr[mc.selected[i]]=='E')
    use_cluster = 0.0;
  else if(secstr[mc.selected[i]]=='L')
    use_cluster = 0.0;
  else if(secstr[mc.selected[i]]=='C')
    use_cluster = USE_CLUSTER;   // 0.1 in define.h
  else
   {
    fprintf(STATUS,"ERROR! secondary structure prediction has a unknown charactrer: %c\n", secstr[mc.selected[i]]);
    exit(1);
   }
  if (mc.selected[i] > nresidues/2.0) {
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1) && native_residue[mc.selected[i]].amino_num!=14) {
	/* can't rotate around a proline phi */
        if (use_cluster > drand48())
         {
          step_phi = desire_phi - cur_phi[mc.selected[i]-1];
          step_phi += GaussianNum()*CLUSTER_NOISE;  // 10. in define.h 
         }
        else
          step_phi = GaussianNum()*CLUSTER_NOISE;

	step_phi *= deg2rad;
//	fprintf(STATUS,"%f\n", step_size);
	
	a = native_residue[mc.selected[i]-1].C;
	b = native_residue[mc.selected[i]].N;
	c = native_residue[mc.selected[i]].CA;
	d = native_residue[mc.selected[i]].C;
        DoRotation(a,b,c,d,step_phi,loop_rotate_natoms[mc.selected[i]][0],loop_rotate_atoms[mc.selected[i]][0]);
        
	mc.delta_phi_angle[i]=step_phi;
        for(j=0;j<loop_rotate_natoms[mc.selected[i]][0];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][0][j]]=1;
        }
      }
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1)) {
        if (use_cluster > drand48())
         {
          step_psi = desire_psi - cur_psi[mc.selected[i]-1];
          step_psi += GaussianNum()*CLUSTER_NOISE;
         }
        else
          step_psi = GaussianNum()*CLUSTER_NOISE;
  
	step_psi *= deg2rad;
	
	a = native_residue[mc.selected[i]].N;
	b = native_residue[mc.selected[i]].CA;
	c = native_residue[mc.selected[i]].C;
	d = native_residue[mc.selected[i]+1].N;
       DoRotation(a,b,c,d,step_psi,loop_rotate_natoms[mc.selected[i]][1],loop_rotate_atoms[mc.selected[i]][1]);
       
	mc.delta_psi_angle[i]=step_psi;
        for(j=0;j<loop_rotate_natoms[mc.selected[i]][1];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][1][j]]=2;
        }
      } 
  }
  else {
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1)) {
        if (use_cluster > drand48())
         {
          step_psi = desire_psi - cur_psi[mc.selected[i]-1];
          step_psi += GaussianNum()*CLUSTER_NOISE;
         }
        else
          step_psi = GaussianNum()*CLUSTER_NOISE;
 
	step_psi *= -deg2rad;
  
	a = native_residue[mc.selected[i]].N;
	c = native_residue[mc.selected[i]].CA;
	b = native_residue[mc.selected[i]].C;
	d = native_residue[mc.selected[i]+1].N;
	DoRotation(a,b,c,d,-step_psi,loop_rotate_natoms[mc.selected[i]][0],loop_rotate_atoms[mc.selected[i]][0]);
  
	mc.delta_psi_angle[i]=-step_psi;
        for(j=0;j<loop_rotate_natoms[mc.selected[i]][0];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][0][j]]=1;
        }
      } 
      if (mc.selected[i] != 0 && mc.selected[i] != (nresidues-1) && native_residue[mc.selected[i]].amino_num!=14) {
        if (use_cluster > drand48())
         {
          step_phi = desire_phi - cur_phi[mc.selected[i]-1];
          step_phi += GaussianNum()*CLUSTER_NOISE;
         }
        else
          step_phi = GaussianNum()*CLUSTER_NOISE;
	
  step_phi *= -deg2rad;
	
	a = native_residue[mc.selected[i]-1].C;
	c = native_residue[mc.selected[i]].N;
	b = native_residue[mc.selected[i]].CA;
	d = native_residue[mc.selected[i]].C;
	DoRotation(a,b,c,d,-step_phi,loop_rotate_natoms[mc.selected[i]][1],loop_rotate_atoms[mc.selected[i]][1]);
  
	mc.delta_phi_angle[i]=-step_phi;
        for(j=0;j<loop_rotate_natoms[mc.selected[i]][1];j++){
          is_rotated[loop_rotate_atoms[mc.selected[i]][1][j]]=2;  
        }
      }
  }
//  check_phipsi();
//   fprintf(STATUS,"%3d %s %3d %3d %7.5f %7.5f %7.5f %7.5f\n",
//       mc.selected[i], native_residue[mc.selected[i]].res, native_residue[mc.selected[i]].amino_num, cluster_bin,
//       desire_phi, cur_phi[mc.selected[i]-1], desire_psi, cur_psi[mc.selected[i]-1]);

  mc.sel_res_num=mc.selected[0];
  UpdateLattice(loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0]);
  all_rotated_natoms = loop_rotate_natoms[mc.selected[i]][0]; // why these two lines?
  all_rotated_atoms = loop_rotate_atoms[mc.selected[i]][0];
  
  NewDeltaContacts(loop_rotate_natoms[mc.selected[i]][0], loop_rotate_atoms[mc.selected[i]][0], loop_not_rotated[mc.selected[i]][0]);
  return;

}
 
void BackboneMove(Float step_size) {
  short temp;
  int a,b,c,d;   /*          d    */
                 /*         /     */ 
                 /*    b - c      */
                 /*   /           */
                 /*  a            */

  mc.delta_angle[0] = step_size;

  if (mc.is_phi) {
    a = native_residue[mc.sel_res_num-1].C;
    b = native_residue[mc.sel_res_num].N;
    c = native_residue[mc.sel_res_num].CA;
    d = native_residue[mc.sel_res_num].C;
    mc.delta_phi_angle[0] = step_size; /* just for proper updating purposes */
    mc.delta_psi_angle[0] = 0;
  }
  else {
    a = native_residue[mc.sel_res_num].N;
    b = native_residue[mc.sel_res_num].CA;
    c = native_residue[mc.sel_res_num].C;
    d = native_residue[mc.sel_res_num+1].N;
    mc.delta_psi_angle[0] = step_size;
    mc.delta_phi_angle[0] = 0;
  } 

  if (mc.sel_res_num <= nresidues/2.0) { 
    mc.delta_angle[0] *=-1;
    mc.delta_psi_angle[0] *=-1;
    mc.delta_phi_angle[0] *=-1;
    temp = b;
    b = c;
    c = temp;
  } 
 
  DoRotation(a,b,c,d,mc.delta_angle[0],rotate_natoms[mc.is_phi][mc.sel_res_num],rotate_atom[mc.is_phi][mc.sel_res_num]);
  //mc.selected[0]=mc.sel_res_num; // added
  all_rotated_natoms = rotate_natoms[mc.is_phi][mc.sel_res_num];
  all_rotated_atoms = rotate_atom[mc.is_phi][mc.sel_res_num];
  
  return;

}

void MakeSidechainMove() {
  int a,b,c,d, i, j;
  int sel_no_chi;
  float p_0to1;
  float cummul_prob;

  mc.sel_rotamer = (int) (drand48()*native_residue[mc.sel_res_num].nrotamers);
  if (USE_ROT_PROB == 1)
   {
    p_0to1 = drand48()*100;
    cummul_prob = 0.;
    for(i=0; i<no_chi_list[native_residue[mc.sel_res_num].amino_num]; i++)
     {
      cummul_prob += prob_ang[native_residue[mc.sel_res_num].amino_num][i]; 
      if (cummul_prob > p_0to1)
	break;
     }
    sel_no_chi = i;
   }
  else
   sel_no_chi = (int) (drand48()*no_chi_list[native_residue[mc.sel_res_num].amino_num]);
  //fprintf(STATUS,"sidechain move at %5d\n", mc.sel_res_num);
  
  old_rotamer = cur_rotamers[mc.sel_res_num];
  cur_rotamers[mc.sel_res_num] = mc.sel_rotamer;

  for(i=0; i<native_residue[mc.sel_res_num].ntorsions; i++) {

    a = sidechain_torsion[mc.sel_res_num][i][0];
    b = sidechain_torsion[mc.sel_res_num][i][1];
    c = sidechain_torsion[mc.sel_res_num][i][2];
    d = sidechain_torsion[mc.sel_res_num][i][3];
//    fprintf(STATUS,"%10ld %3d %s chi: %d  no_list: %2d %8.3f %8.3f %8.3f %8.3f %8.3f\n", mcstep, mc.sel_res_num, native_residue[mc.sel_res_num].res, i,
//           sel_no_chi, rotamer_angles[mc.sel_res_num].chis[sel_no_chi][i], deviation_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi][i],
//	   p_0to1, cummul_prob, prob_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi]);

    if (USE_ROTAMERS)
      mc.delta_angle[i] = rotamer_angles[mc.sel_res_num].chis[sel_no_chi][i] +
                     	  GaussianNum()*deviation_ang[native_residue[mc.sel_res_num].amino_num][sel_no_chi][i]
                  	  - native_residue[mc.sel_res_num].chi[i]; 
    else
      mc.delta_angle[i] = GaussianNum()*SIDECHAIN_NOISE;
    native_residue[mc.sel_res_num].tmpchi[i] = native_residue[mc.sel_res_num].chi[i] + mc.delta_angle[i];
    
    DoRotation(a,b,c,d,mc.delta_angle[i],rotate_sidechain_natoms[mc.sel_res_num][i],rotate_sidechain_atom[mc.sel_res_num][i]);
    for(j=0;j<rotate_sidechain_natoms[mc.sel_res_num][i];j++){
      is_rotated[rotate_sidechain_atom[mc.sel_res_num][i][j]]=i+1;
    }
  }

  all_rotated_natoms = rotate_sidechain_natoms[mc.sel_res_num][0]; 
  all_rotated_atoms = rotate_sidechain_atom[mc.sel_res_num][0];
  
  return;
}

void SidechainMove() {

  all_rotated_natoms = 0;
  total_pairs=total_pairs2=0;  
  total_hbond_pairs=0;

  do {
    mc.sel_res_num = (int) (drand48()*nresidues);
  } while (native_residue[mc.sel_res_num].nrotamers <=1 || native_residue[mc.sel_res_num].amino_num==14);

  MakeSidechainMove();
  sidemovedone = 1;
  UpdateLattice(rotate_sidechain_natoms[mc.sel_res_num][0],rotate_sidechain_atom[mc.sel_res_num][0]);   
  NewDeltaContacts(rotate_sidechain_natoms[mc.sel_res_num][0],rotate_sidechain_atom[mc.sel_res_num][0],sidechain_not_rotated[mc.sel_res_num][0]);
  
  return;
}


