void InitializeFraglib();
void InitializeFragRotationData();
Float RandomAngleGenerate(Float miu, Float sigma);
void FragMove();
void CrankFragMove();


void InitializeFraglib()
{
	FILE *ffraglib;
  	char line[200];
  	float phi,psi;
  	int i,j,k;
  	int posi_cnt,cand_cnt,cand_len;
  	posi_cnt=cand_cnt=cand_len=0;
  	
  	fraglib_ncands=(int *) calloc(nresidues-6, sizeof(int));
  	fraglib_cands=(float ***) calloc(nresidues-6, sizeof(float **));
  	
  	for (i=0;i<nresidues-6;i++)
  	{
  		fraglib_cands[i]=(float **) calloc(1000, sizeof(float *));  // max:1000 candidates in each position
  		for (j=0;j<1000;j++)
  		{
  			fraglib_cands[i][j]=(float *) calloc(21,sizeof(float));  // the first one is length; then up to 10*2 angle values
  		}
  	}
  	
  	fprintf(STATUS,"Opening the file: %s\n", fraglib_file);
  	
  	if((ffraglib=fopen(fraglib_file,"r"))==NULL)
   	{
    	fprintf(STATUS,"ERROR: Can't open the file: %s!\n", fraglib_file);
    	exit(1);
   	}
   	
  	while(fgets(line,150,ffraglib)!=NULL)
  	{
  		if (!strncmp(line,"position",8))
  		{
  			cand_len=0;
  			cand_cnt=-1;
  			sscanf(line,"%*s %d %*s %d",&i,&j);
  			posi_cnt=i;
  			if (j<=1000)
  				fraglib_ncands[i]=j;
  			else
  				fraglib_ncands[i]=1000;
  		}
  		else if (line[0]!='\n' && line[0]!='\r')
  		{
  			if (cand_cnt<1000)
  			{
  				sscanf(line,"%*s %*s %*d %*s %*s %f %f %*f",&phi,&psi);
  				fraglib_cands[posi_cnt][cand_cnt][cand_len*2+1]=phi;
  				fraglib_cands[posi_cnt][cand_cnt][cand_len*2+1+1]=psi;
  				cand_len++;
  			}
  		}
  		else
  		{	
  			if (cand_cnt<1000)
  			{
  				if (cand_cnt!=-1)
  					fraglib_cands[posi_cnt][cand_cnt][0]=(float)cand_len;
  				cand_len=0;
  				cand_cnt++;
  			}
  		}
  	}
  	
  	/* for test 
	for (i=0;i<nresidues-6;i++)
	{
		printf ("Posi:%d  Num:%d\n",i,(int)fraglib_ncands[i]);
		for (j=0;j<fraglib_ncands[i];j++)
		{
			printf ("\nLength:%d\n",(int)fraglib_cands[i][j][0]);
			for (k=1;k<fraglib_cands[i][j][0]*2+1;k+=2)
				printf ("%8.3f %8.3f\n",fraglib_cands[i][j][k],fraglib_cands[i][j][k+1]);
		}
	}
	printf("Read Fragment Library Successfully\n");
	end test */
  	
  	fprintf(STATUS,"Read Fragment Library Successfully: %s\n",fraglib_file);
  	return; 		
}



void InitializeFragRotationData()
{
	int i,j;
	// rotate_atom[0=psi, 1=phi][which residue][list of atoms] 
	frag_left_rotate_natoms = (short **) calloc(2, sizeof(short *));     //2*n1    where n1:nresidues n2:natoms 
  	frag_left_rotate_atom = (short ***) calloc(2, sizeof(short **));     // 2*n1*n2
  	frag_left_not_rotated = (char ***) calloc(2,sizeof(char **));        // 2*n1*n2

  	frag_right_rotate_natoms = (short **) calloc(2, sizeof(short *));     //2*n1    where n1:nresidues n2:natoms 
  	frag_right_rotate_atom = (short ***) calloc(2, sizeof(short **));     // 2*n1*n2
  	frag_right_not_rotated = (char ***) calloc(2,sizeof(char **));        // 2*n1*n2


  	for(i=0; i<2; i++) 
  	{
    	frag_left_rotate_atom[i] = (short **) calloc(nresidues, sizeof(short *));
    	frag_left_not_rotated[i] = (char **) calloc(nresidues, sizeof(char *));
    	frag_left_rotate_natoms[i] = (short *) calloc(nresidues,sizeof(short));

    	frag_right_rotate_atom[i] = (short **) calloc(nresidues, sizeof(short *));
    	frag_right_not_rotated[i] = (char **) calloc(nresidues, sizeof(char *));
    	frag_right_rotate_natoms[i] = (short *) calloc(nresidues,sizeof(short));
    
    	for(j=0; j<nresidues; j++) 
    	{
      		frag_left_rotate_atom[i][j] = (short *) calloc(natoms, sizeof(short));
      		frag_left_not_rotated[i][j] = (char *) calloc(natoms,sizeof(char));
    	}

    	for(j=0; j<nresidues; j++) 
    	{
      		frag_right_rotate_atom[i][j] = (short *) calloc(natoms, sizeof(short));
      		frag_right_not_rotated[i][j] = (char *) calloc(natoms,sizeof(char));
    	}
  	}


  	for(i=0; i<nresidues; i++)   // like i<=nresidues/2 left rotate
  	{
  		frag_left_rotate_natoms[0][i]=0;
    	frag_left_rotate_natoms[1][i]=0;

    	for(j=0; j<natoms; j++)
      	{   
			if (native[j].res_num < i) 
			{
	  			frag_left_rotate_atom[0][i][frag_left_rotate_natoms[0][i]++] = j;
	  			frag_left_rotate_atom[1][i][frag_left_rotate_natoms[1][i]++] = j;
			}
			else if (native[j].res_num == i) 
			{
	  			if (j != native_residue[i].O && j != native_residue[i].C && j != native_residue[i].CA)
	    			frag_left_rotate_atom[0][i][frag_left_rotate_natoms[0][i]++] = j;  // psi
			}
		}
  	}

  	for(i=0; i<nresidues; i++)   // like i>nresidues/2 right rotate
  	{
    	frag_right_rotate_natoms[0][i]=0;
    	frag_right_rotate_natoms[1][i]=0;
    
      	for(j=0; j<natoms; j++)
      	{
      		if (native[j].res_num > i) 
      		{
      			frag_right_rotate_atom[0][i][frag_right_rotate_natoms[0][i]++] = j; 
	  			frag_right_rotate_atom[1][i][frag_right_rotate_natoms[1][i]++] = j; 
      		}	
			else if (native[j].res_num == i) 
			{
	  			if (j != native_residue[i].N && j!= native_residue[i].CA)
	    			frag_right_rotate_atom[1][i][frag_right_rotate_natoms[1][i]++] = j;  // phi
	  			if (j == native_residue[i].O)
	    			frag_right_rotate_atom[0][i][frag_right_rotate_natoms[0][i]++] = j;  // psi
			}
		}
	}

	       
   
  	for(i=0; i<nresidues; i++) 
  	{
    	for(j=0; j<frag_left_rotate_natoms[0][i]; j++)
      		frag_left_not_rotated[0][i][frag_left_rotate_atom[0][i][j]]=1;
    	for(j=0; j<frag_left_rotate_natoms[1][i]; j++)
      		frag_left_not_rotated[1][i][frag_left_rotate_atom[1][i][j]]=1;

      	for(j=0; j<frag_right_rotate_natoms[0][i]; j++)
      		frag_right_not_rotated[0][i][frag_right_rotate_atom[0][i][j]]=1;
    	for(j=0; j<frag_right_rotate_natoms[1][i]; j++)
      		frag_right_not_rotated[1][i][frag_right_rotate_atom[1][i][j]]=1;

    	for(j=0; j<natoms; j++) 
    	{
      		frag_left_not_rotated[0][i][j]=!frag_left_not_rotated[0][i][j];
      		frag_left_not_rotated[1][i][j]=!frag_left_not_rotated[1][i][j];            //   0:rotated 1:not rotate
      		frag_right_not_rotated[0][i][j]=!frag_right_not_rotated[0][i][j];
      		frag_right_not_rotated[1][i][j]=!frag_right_not_rotated[1][i][j];
    	}
  	}
}



Float RandomAngleGenerate(Float miu, Float sigma)
{
	return GaussianNum()*sigma+miu;
}



void FragMove()
{
	int sel_idx,fraglength;
	int a,b,c,d,i,j,k;
	double rotate_phi,rotate_psi,desire_phi,desire_psi;
	//double bef_phi[MAXSEQUENCE], bef_psi[MAXSEQUENCE];

	
  	all_rotated_natoms = 0;
  	total_pairs=total_pairs2=0;
  	total_hbond_pairs=0;
  	nomove = 0;

  	mc.sel_res_num = (int)(drand48()*(nresidues-6));  // minlength=7
  	sel_idx=(int)(drand48()*fraglib_ncands[mc.sel_res_num]);
  	fraglength=fraglib_cands[mc.sel_res_num][sel_idx][0];  	  	
  	
  	mc.loop_size = fraglength;


  	for (i=mc.sel_res_num;i<mc.sel_res_num+fraglength;i++)
  	{
  		if(is_template[i]==1)
   		{
    		nomove = 1;
    		nothers++;
    		break;
   		}
  	}
  	if (nomove==1)
  	{
  		return;
  	}

  	/*test
  	check_phipsi();

  	int tmpi;
  	for (tmpi=0;tmpi<nresidues-2;tmpi++)
    {
        bef_phi[tmpi]=cur_phi[tmpi];
        bef_psi[tmpi]=cur_psi[tmpi];
    }
    //end of test*/



  	if (mc.sel_res_num+fraglength-1>nresidues/2)
  	{
  		for (i=mc.sel_res_num;i<mc.sel_res_num+fraglength;i++)
  		{
      		// phi
      		if (i!=0 && native_residue[i].amino_num!=14)
      		{
      			desire_phi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+1];
      			//rotate_phi=(desire_phi-cur_phi[i-1])*deg2rad;
      			rotate_phi=(RandomAngleGenerate(desire_phi,NOISE_RANGE_PHI)-cur_phi[i-1])*deg2rad;

      			a = native_residue[i-1].C;
      			b = native_residue[i].N;
				c = native_residue[i].CA;
				d = native_residue[i].C;
				
        		DoRotation(a,b,c,d,rotate_phi,frag_right_rotate_natoms[1][i],frag_right_rotate_atom[1][i]);
        		
        		mc.delta_phi_angle[i-mc.sel_res_num]=rotate_phi;
        		
        		for(j=0;j<frag_right_rotate_natoms[1][i];j++)
        		{
          			is_rotated[frag_right_rotate_atom[1][i][j]]+=1;
        		}
      		}
      		
      		//psi
      		if (i!=nresidues-1)
      		{
      			desire_psi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+2];
      			//rotate_psi=(desire_psi-cur_psi[i-1])*deg2rad;
      			rotate_psi=(RandomAngleGenerate(desire_psi,NOISE_RANGE_PSI)-cur_psi[i-1])*deg2rad;
      			a = native_residue[i].N;
      			b = native_residue[i].CA;
				c = native_residue[i].C;
				d = native_residue[i+1].N;
				
        		DoRotation(a,b,c,d,rotate_psi,frag_right_rotate_natoms[0][i],frag_right_rotate_atom[0][i]);
        		
        		mc.delta_psi_angle[i-mc.sel_res_num]=rotate_psi;
        		
        		for(j=0;j<frag_right_rotate_natoms[0][i];j++)
        		{
          			is_rotated[frag_right_rotate_atom[0][i][j]]+=1;
        		}
      		}
  		}
  	}

  	else 
  	{   
  		for (i=mc.sel_res_num+fraglength-1;i>mc.sel_res_num-1;i--)
  		{
  			//psi
  			if (i!=nresidues-1)
      		{
      			desire_psi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+2];
      			//rotate_psi=(desire_psi-cur_psi[i-1])*(-deg2rad);
     			rotate_psi=(RandomAngleGenerate(desire_psi,NOISE_RANGE_PSI)-cur_psi[i-1])*(-deg2rad);

      			a = native_residue[i].N;
      			b = native_residue[i].CA;
				c = native_residue[i].C;
				d = native_residue[i+1].N;
        		DoRotation(a,b,c,d,rotate_psi,frag_left_rotate_natoms[0][i],frag_left_rotate_atom[0][i]);

        		mc.delta_psi_angle[i-mc.sel_res_num]=rotate_psi;
        		for(j=0;j<frag_left_rotate_natoms[0][i];j++)
        		{
          			is_rotated[frag_left_rotate_atom[0][i][j]]+=1;
        		}
      		}

      		//phi
      		if (i!=0 && native_residue[i].amino_num!=14)
      		{
      			desire_phi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+1];
      			//rotate_phi=(desire_phi-cur_phi[i-1])*(-deg2rad);
      			rotate_phi=(RandomAngleGenerate(desire_phi,NOISE_RANGE_PHI)-cur_phi[i-1])*(-deg2rad);
      			
      			a = native_residue[i-1].C;
      			b = native_residue[i].N;
				c = native_residue[i].CA;
				d = native_residue[i].C;
        		DoRotation(a,b,c,d,rotate_phi,frag_left_rotate_natoms[1][i],frag_left_rotate_atom[1][i]);

        		mc.delta_phi_angle[i-mc.sel_res_num]=rotate_phi;
        		for(j=0;j<frag_left_rotate_natoms[1][i];j++)
        		{
          			is_rotated[frag_left_rotate_atom[1][i][j]]+=1;
        		}
      		} 		
  		}
  	}

  	/*test
  	
    check_phipsi();
    printf ("Phi and Psi before/after the move\n");
    for(tmpi=0;tmpi<nresidues-2;tmpi++)
    {
        if (abs(bef_phi[tmpi]-cur_phi[tmpi])==0 && abs(bef_psi[tmpi]-cur_psi[tmpi])==0)
            ;//printf ("phi: %8.3f %8.3f  psi: %8.3f %8.3f\n",bef_phi[tmpi],cur_phi[tmpi],bef_psi[tmpi],cur_psi[tmpi]);
        else
            printf ("phi: %8.3f %8.3f  psi: %8.3f %8.3f *****  res_idx:%d res_sel_idx:%d desire_phi:%8.3f desire_psi:%8.3f\n"\
            	,bef_phi[tmpi],cur_phi[tmpi],bef_psi[tmpi],cur_psi[tmpi],tmpi,mc.sel_res_num,\
            	fraglib_cands[mc.sel_res_num][sel_idx][(tmpi+1-mc.sel_res_num)*2+1],fraglib_cands[mc.sel_res_num][sel_idx][(tmpi+1-mc.sel_res_num)*2+2]);
    } 
    //printf ("A Fragment Move Finish\n");//end of test*/


  	

  	if (mc.sel_res_num+fraglength-1>nresidues/2)
  	{
  		all_rotated_natoms = frag_right_rotate_natoms[1][mc.sel_res_num]; 
  		all_rotated_atoms = frag_right_rotate_atom[1][mc.sel_res_num];
  		UpdateLattice(frag_right_rotate_natoms[1][mc.sel_res_num], frag_right_rotate_atom[1][mc.sel_res_num]);	
  		NewDeltaContacts(frag_right_rotate_natoms[1][mc.sel_res_num], frag_right_rotate_atom[1][mc.sel_res_num], frag_right_not_rotated[1][mc.sel_res_num]);
  		return;
  	}
  	else
  	{
  		all_rotated_natoms = frag_left_rotate_natoms[0][mc.sel_res_num+fraglength-1]; 
  		all_rotated_atoms = frag_left_rotate_atom[0][mc.sel_res_num+fraglength-1];
  		UpdateLattice(frag_left_rotate_natoms[0][mc.sel_res_num+fraglength-1], frag_left_rotate_atom[0][mc.sel_res_num+fraglength-1]); 	
  		NewDeltaContacts(frag_left_rotate_natoms[0][mc.sel_res_num+fraglength-1], frag_left_rotate_atom[0][mc.sel_res_num+fraglength-1], frag_left_not_rotated[0][mc.sel_res_num+fraglength-1]);
  		return;
  	}
}



void CrankFragMove()
{
	int sel_idx,fraglength;
	int a,b,c,d,i,j,k;
	double rotate_phi,rotate_psi,desire_phi,desire_psi;
	//double bef_phi[MAXSEQUENCE], bef_psi[MAXSEQUENCE];

	
  	all_rotated_natoms = 0;
  	total_pairs=total_pairs2=0;
  	total_hbond_pairs=0;
  	nomove = 0;

  	mc.sel_res_num = (int)(drand48()*(nresidues-6));  // minlength=7
  	sel_idx=(int)(drand48()*fraglib_ncands[mc.sel_res_num]);
  	fraglength=fraglib_cands[mc.sel_res_num][sel_idx][0];  	  	
  	
  	mc.loop_size = fraglength;


  	for (i=mc.sel_res_num;i<mc.sel_res_num+fraglength;i++)
  	{
  		if(is_template[i]==1)
   		{
    		nomove = 1;
    		nothers++;
    		break;
   		}
  	}
  	if (nomove==1)
  	{
  		return;
  	}

  	/*test
  	check_phipsi();

  	int tmpi;
  	for (tmpi=0;tmpi<nresidues-2;tmpi++)
    {
        bef_phi[tmpi]=cur_phi[tmpi];
        bef_psi[tmpi]=cur_psi[tmpi];
    }
    //end of test*/



  	if (mc.sel_res_num+fraglength-1>nresidues/2)
  	{
  		for (i=mc.sel_res_num;i<mc.sel_res_num+fraglength;i++)
  		{
      		// phi
      		if (i!=0 && native_residue[i].amino_num!=14)
      		{
      			if (i!=mc.sel_res_num) // omits the first phi
      			{
      				desire_phi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+1];
	      			rotate_phi=(RandomAngleGenerate(desire_phi,NOISE_RANGE_PHI)-cur_phi[i-1])*deg2rad;

	      			a = native_residue[i-1].C;
	      			b = native_residue[i].N;
					c = native_residue[i].CA;
					d = native_residue[i].C;
					
	        		DoRotation(a,b,c,d,rotate_phi,frag_right_rotate_natoms[1][i],frag_right_rotate_atom[1][i]);
	        		
	        		mc.delta_phi_angle[i-mc.sel_res_num]=rotate_phi;
	        		
	        		for(j=0;j<frag_right_rotate_natoms[1][i];j++)
	        		{
	          			is_rotated[frag_right_rotate_atom[1][i][j]]+=1;
	        		}
      			}	
      		}
      		
      		//psi
      		if (i!=nresidues-1)
      		{
      			if (i!=mc.sel_res_num+fraglength-1) // omits the last psi
      			{
      				desire_psi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+2];
	      			rotate_psi=(RandomAngleGenerate(desire_psi,NOISE_RANGE_PSI)-cur_psi[i-1])*deg2rad;
	      			a = native_residue[i].N;
	      			b = native_residue[i].CA;
					c = native_residue[i].C;
					d = native_residue[i+1].N;
					
	        		DoRotation(a,b,c,d,rotate_psi,frag_right_rotate_natoms[0][i],frag_right_rotate_atom[0][i]);
	        		
	        		mc.delta_psi_angle[i-mc.sel_res_num]=rotate_psi;
	        		
	        		for(j=0;j<frag_right_rotate_natoms[0][i];j++)
	        		{
	          			is_rotated[frag_right_rotate_atom[0][i][j]]+=1;
	        		}
      			}
      		}
  		}
  	}

  	else 
  	{   
  		for (i=mc.sel_res_num+fraglength-1;i>mc.sel_res_num-1;i--)
  		{
  			//psi
  			if (i!=nresidues-1)
      		{
      			if (i!=mc.sel_res_num+fraglength-1)
      			{
      				desire_psi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+2];
	     			rotate_psi=(RandomAngleGenerate(desire_psi,NOISE_RANGE_PSI)-cur_psi[i-1])*(-deg2rad);

	      			a = native_residue[i].N;
	      			b = native_residue[i].CA;
					c = native_residue[i].C;
					d = native_residue[i+1].N;
	        		DoRotation(a,b,c,d,rotate_psi,frag_left_rotate_natoms[0][i],frag_left_rotate_atom[0][i]);

	        		mc.delta_psi_angle[i-mc.sel_res_num]=rotate_psi;
	        		for(j=0;j<frag_left_rotate_natoms[0][i];j++)
	        		{
	          			is_rotated[frag_left_rotate_atom[0][i][j]]+=1;
	        		}
      			}
      		}

      		//phi
      		if (i!=0 && native_residue[i].amino_num!=14)
      		{
      			if (i!=mc.sel_res_num)
      			{
      				desire_phi=fraglib_cands[mc.sel_res_num][sel_idx][(i-mc.sel_res_num)*2+1];
	      			rotate_phi=(RandomAngleGenerate(desire_phi,NOISE_RANGE_PHI)-cur_phi[i-1])*(-deg2rad);
	      			
	      			a = native_residue[i-1].C;
	      			b = native_residue[i].N;
					c = native_residue[i].CA;
					d = native_residue[i].C;
	        		DoRotation(a,b,c,d,rotate_phi,frag_left_rotate_natoms[1][i],frag_left_rotate_atom[1][i]);

	        		mc.delta_phi_angle[i-mc.sel_res_num]=rotate_phi;
	        		for(j=0;j<frag_left_rotate_natoms[1][i];j++)
	        		{
	          			is_rotated[frag_left_rotate_atom[1][i][j]]+=1;
	        		}
      			}
      		} 		
  		}
  	}

  	/*test
  	
    check_phipsi();
    printf ("Phi and Psi before/after the move\n");
    for(tmpi=0;tmpi<nresidues-2;tmpi++)
    {
        if (abs(bef_phi[tmpi]-cur_phi[tmpi])==0 && abs(bef_psi[tmpi]-cur_psi[tmpi])==0)
            ;//printf ("phi: %8.3f %8.3f  psi: %8.3f %8.3f\n",bef_phi[tmpi],cur_phi[tmpi],bef_psi[tmpi],cur_psi[tmpi]);
        else
            printf ("phi: %8.3f %8.3f  psi: %8.3f %8.3f *****  res_idx:%d res_sel_idx:%d desire_phi:%8.3f desire_psi:%8.3f\n"\
            	,bef_phi[tmpi],cur_phi[tmpi],bef_psi[tmpi],cur_psi[tmpi],tmpi,mc.sel_res_num,\
            	fraglib_cands[mc.sel_res_num][sel_idx][(tmpi+1-mc.sel_res_num)*2+1],fraglib_cands[mc.sel_res_num][sel_idx][(tmpi+1-mc.sel_res_num)*2+2]);
    } 
    //printf ("A Fragment Move Finish\n");//end of test*/
  	

  	if (mc.sel_res_num+fraglength-1>nresidues/2)
  	{
  		all_rotated_natoms = frag_right_rotate_natoms[0][mc.sel_res_num]; 
  		all_rotated_atoms = frag_right_rotate_atom[0][mc.sel_res_num];
  		UpdateLattice(frag_right_rotate_natoms[0][mc.sel_res_num], frag_right_rotate_atom[0][mc.sel_res_num]);	
  		NewDeltaContacts(frag_right_rotate_natoms[0][mc.sel_res_num], frag_right_rotate_atom[0][mc.sel_res_num], frag_right_not_rotated[0][mc.sel_res_num]);
  		return;
  	}
  	else
  	{
  		all_rotated_natoms = frag_left_rotate_natoms[1][mc.sel_res_num+fraglength-1]; 
  		all_rotated_atoms = frag_left_rotate_atom[1][mc.sel_res_num+fraglength-1];
  		UpdateLattice(frag_left_rotate_natoms[1][mc.sel_res_num+fraglength-1], frag_left_rotate_atom[1][mc.sel_res_num+fraglength-1]); 	
  		NewDeltaContacts(frag_left_rotate_natoms[1][mc.sel_res_num+fraglength-1], frag_left_rotate_atom[1][mc.sel_res_num+fraglength-1], frag_left_not_rotated[1][mc.sel_res_num+fraglength-1]);
  		return;
  	}
}








