/*
 *  Utils.cpp
 *  Created by Jeremy Riousset on 11/19/07.
 */

#include "MpiFunctions.h"

void MPI_foo::CreateComm(void)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Var::world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_Var::n_processes);
};

void MPI_foo::CreateGridComm(void)
{
	//	int processes_per_dimensions = (int)pow(MPI_Var::n_processes,1/3.);
	MPI_Var::dim_sizes		= MPI_Var::n_processes;										//processes_per_dimensions;
	MPI_Var::wrap_around	= 0;														// The domain is NOT periodic in z
	int tmp_grid_rank;
	
	MPI_Cart_create(MPI_COMM_WORLD, MPI_Var::n_dims, &MPI_Var::dim_sizes, &MPI_Var::wrap_around, MPI_Var::reorder, &MPI_Var::grid_comm);
	MPI_Comm_rank(MPI_Var::grid_comm, &tmp_grid_rank);
	MPI_Cart_coords(MPI_Var::grid_comm, tmp_grid_rank, MPI_Var::max_dims, &MPI_Var::coordinates);
	MPI_Cart_rank(MPI_Var::grid_comm, &MPI_Var::coordinates, &MPI_Var::grid_rank);
};

void MPI_foo::CreateCartComm(void)
{
	MPI_Var::free_coords = 1;
    MPI_Cart_sub(	MPI_Var::grid_comm, &MPI_Var::free_coords, &MPI_Var::x_comm);
	MPI_Comm_rank(	MPI_Var::x_comm,	&MPI_Var::x_rank);
	MPI_Cart_shift(	MPI_Var::x_comm,	MPI_Var::direction, MPI_Var::displacement, &MPI_Var::prev_x_rank, &MPI_Var::next_x_rank);
};

void MPI_foo::InitLocalDimensions(void)
{
	Var::local_N.x = (Var::N.x-2)/MPI_Var::dim_sizes+2;
	Var::local_N.y = Var::N.y;
	Var::local_N.z = Var::N.z;
};

void MPI_foo::FreeComm(void)
{
	MPI_Comm_free( &MPI_Var::grid_comm);
	MPI_Comm_free( &MPI_Var::x_comm);
};

CMatrix3D MPI_foo::Scatter(CMatrix3D MM)
{
	CMatrix3D	llocal_MM;
	int			source;
	int			dest;
	int			tag		=0;
	MPI_Status	status;
	
	// SCATTER 1-D //
	llocal_MM.init(Var::local_N.x,Var::N.y,Var::N.z);
	if(MPI_Var::x_rank == MPI_Var::root)
	{
		//COPY global plane MM(0,jj,kk)     to tmp_M1(0,jj,kk)	   on first process through x_comm//			
		memcpy(&llocal_MM.pElems[0], &MM.pElems[0], Var::N.y*Var::N.z*sizeof(double));
		//SEND global plane MM(N.x-1,jj,kk) to tmp_M1(N.x-1,jj,kk) on last  process through x_comm//
		dest = MPI_Var::dim_sizes-1;
		MPI_Send(	&MM.pElems[(Var::N.x-1)*Var::N.y*Var::N.z], Var::N.y*Var::N.z,	MPI_DOUBLE, dest, tag, MPI_Var::x_comm);
	}
	if(MPI_Var::x_rank == MPI_Var::dim_sizes-1)
	{
		//RECEIVE global plane MM(N.x-1,jj,kk) to tmp_M1(N.x-1,jj,kk) on last  process through x_comm//
		source = MPI_Var::root;
		MPI_Recv(	&llocal_MM(Var::local_N.x-1,0,0), Var::N.y*Var::N.z, MPI_DOUBLE, source, tag, MPI_Var::x_comm, &status);
	}
	//SCATTER all intermediate plane through x_comm//
	MPI_Scatter(&MM.pElems[Var::N.y*Var::N.z], (Var::local_N.x-2)*Var::N.y*Var::N.z, MPI_DOUBLE, 
				&llocal_MM.pElems[Var::N.y*Var::N.z], (Var::local_N.x-2)*Var::N.y*Var::N.z, MPI_DOUBLE, 
				MPI_Var::root, MPI_Var::x_comm);
	// END OF SCATTER 1-D //
	return llocal_MM;
};

Charge MPI_foo::Scatter(Charge CC)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CHARGE TYPE TO SAVE COMMUNICATION (BCAST) TIME //
	Charge local_CC;
	local_CC.rho  = Scatter(CC.rho);
	local_CC.Un   = Scatter(CC.Un);
	return local_CC;
};

Potential MPI_foo::Scatter(Potential PP)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CHARGE TYPE TO SAVE COMMUNICATION (BCAST) TIME //
	Potential local_PP;
	if (MPI_Var::world_rank == MPI_Var::root)
	{
		local_PP.EquiPotential	= PP.EquiPotential;
		local_PP.Vo				= PP.Vo;
	}
	MPI_Bcast(&local_PP.EquiPotential, 1, MPI_INT,		0, MPI_COMM_WORLD);
	MPI_Bcast(&local_PP.Vo,			1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		
	local_PP.rho	= Scatter(PP.rho);
	local_PP.Un		= Scatter(PP.Un);
	return local_PP;
};

CMatrix3D MPI_foo::Gather(CMatrix3D llocal_MM)
{
	CMatrix3D MM;
	int			source;
	int			dest;
	int			tag		= 0;
	MPI_Status	status;
	
	if(MPI_Var::x_rank == MPI_Var::root) MM.init(Var::N.x,Var::N.y,Var::N.z);
	if(MPI_Var::x_rank == MPI_Var::dim_sizes-1)
	{
		//RECEIVE global plane MM(N.x-1,jj,kk) to tmp_M1(N.x-1,jj,kk) on last  process through x_comm//
		dest = MPI_Var::root;
		MPI_Send(	&llocal_MM(Var::local_N.x-1,0,0), Var::N.y*Var::N.z, MPI_DOUBLE, dest, tag, MPI_Var::x_comm);
	}
	if(MPI_Var::x_rank == MPI_Var::root)
	{
		//COPY global plane MM(0,jj,kk)     to tmp_M1(0,jj,kk)	   on first process through x_comm//			
		memcpy(&MM.pElems[0], &llocal_MM.pElems[0], Var::N.y*Var::N.z*sizeof(double));
		//SEND global plane MM(N.x-1,jj,kk) to tmp_M1(N.x-1,jj,kk) on last  process through x_comm//
		source = MPI_Var::dim_sizes-1;
		MPI_Recv(	&MM.pElems[(Var::N.x-1)*Var::N.y*Var::N.z], Var::N.y*Var::N.z,	MPI_DOUBLE, source, tag, MPI_Var::x_comm, &status);
	}
	//SCATTER all intermediate plane through x_comm//
	MPI_Gather(&llocal_MM.pElems[Var::N.y*Var::N.z], (Var::local_N.x-2)*Var::N.y*Var::N.z, MPI_DOUBLE, 
			   &MM.pElems[Var::N.y*Var::N.z], (Var::local_N.x-2)*Var::N.y*Var::N.z, MPI_DOUBLE, 
			   MPI_Var::root, MPI_Var::x_comm);
	return MM;
};

void MPI_foo::CreateGhostPlane(void)
{
	int				nb_elems_in_x_plane = (Var::local_N.y-2)*(Var::local_N.z-2)/2;
	int				displacements_red[nb_elems_in_x_plane];
	int				displacements_blk[nb_elems_in_x_plane];
	int				block_lengths_red[nb_elems_in_x_plane];
	int				block_lengths_blk[nb_elems_in_x_plane];
	int				n_red, n_blk;

	MPI_Var::front_plane.init(Var::N.y,Var::N.z);
	MPI_Var::back_plane.init(Var::N.y,Var::N.z);
	MPI_Var::left_plane.init(Var::N.x-2,Var::N.z);
	MPI_Var::right_plane.init(Var::N.x-2,Var::N.z);
	MPI_Var::top_plane.init(Var::N.x-2,Var::N.y-2);
	MPI_Var::bot_plane.init(Var::N.x-2,Var::N.y-2);

	MPI_Type_contiguous(Var::N.y*Var::N.z, MPI_DOUBLE, &MPI_Var::x_plane);
	MPI_Type_commit(&MPI_Var::x_plane);
	MPI_Type_vector(Var::N.x-2, Var::N.z, Var::N.y*Var::N.z, MPI_DOUBLE, &MPI_Var::y_plane);
	MPI_Type_commit(&MPI_Var::y_plane);
	MPI_Type_vector((Var::N.x-2)*(Var::N.y-2), 1, Var::N.z, MPI_DOUBLE, &MPI_Var::z_plane);
	MPI_Type_commit(&MPI_Var::z_plane);

	/* CREATE LOCAL RED AND BLACK X-PLANES */
	n_red = 0;
	n_blk = 0;
	for (int jj = 1 ; jj <Var::local_N.y-1 ; jj++)	for (int kk = 1 ; kk <Var::local_N.z-1 ; kk++)
	{
		if( (jj+kk)%2 == 0)	
		{	
			displacements_blk[n_blk] = jj*Var::local_N.z + kk;
			block_lengths_blk[n_blk] = 1;
			n_blk++;
		}
		if( (jj+kk)%2 == 1)	
		{	
			displacements_red[n_red] = jj*Var::local_N.z + kk;
			block_lengths_red[n_red] = 1;
			n_red++;
		}
	}
	MPI_Type_indexed(nb_elems_in_x_plane, block_lengths_red, displacements_red, MPI_DOUBLE,&MPI_Var::local_x_plane_red);
	MPI_Type_indexed(nb_elems_in_x_plane, block_lengths_blk, displacements_blk, MPI_DOUBLE,&MPI_Var::local_x_plane_blk);
	MPI_Type_commit(&MPI_Var::local_x_plane_red);	
	MPI_Type_commit(&MPI_Var::local_x_plane_blk);	
};

void	MPI_foo::Reduce(CMatrix3D&	MM)
{
	MPI_Reduce(&MPI_Var::back_plane[0][0],	&MM.pElems[0],								Var::N.y*Var::N.z, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
	MPI_Reduce(&MPI_Var::front_plane[0][0],	&MM.pElems[(Var::N.x-1)*Var::N.y*Var::N.z],	Var::N.y*Var::N.z, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
	for(int ii = 1 ; ii<Var::N.x-1 ; ii++)
	{
		MPI_Reduce(&MPI_Var::left_plane[ii-1][0],	&MM.pElems[ii*Var::N.y*Var::N.z],							Var::N.z, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
		MPI_Reduce(&MPI_Var::right_plane[ii-1][0],	&MM.pElems[ii*Var::N.y*Var::N.z+(Var::N.y-1)*Var::N.z],		Var::N.z, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
		for(int jj = 1 ; jj<Var::N.y-1 ; jj++)
		{
			MPI_Reduce(&MPI_Var::bot_plane[ii-1][jj-1],	&MM.pElems[ii*Var::N.y*Var::N.z+jj*Var::N.z],				1, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
			MPI_Reduce(&MPI_Var::top_plane[ii-1][jj-1],	&MM.pElems[ii*Var::N.y*Var::N.z+jj*Var::N.z+Var::N.z-1],	1, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::x_comm);
		}
	}
};
