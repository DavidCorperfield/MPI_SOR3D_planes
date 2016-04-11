/*
 *  MPI_Input.h
 *  Created by Jeremy Riousset on 11/19/07.
 */

#ifndef MPIINPUT_H
#define MPIINPUT_H

#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#include "Matrix.h"

using namespace std;

/*MPI VARIABLES DECLARATION********************************************************/
class MPI_Var
{
public:
	/*CREATE COMM*/
	static	int						root;
	static	int						world_rank;
	static	int						n_processes;
	
	/*CREATE GRID COMMUNICATOR*/
	static	int						grid_rank;
	static	int						n_dims;
	static	int						dim_sizes;
	static	int						max_dims;
	static	int						coordinates;
	static	int						wrap_around;
	static	int						reorder;
	static	MPI_Comm				grid_comm;
	
	/*CREATE CARTESIAN COMMUNICATOR*/
	static	int						free_coords;
	static	int						direction;
	static	int						displacement;
	static	int						prev_x_rank;
	static	int							 x_rank;
	static	int						next_x_rank;
	static	MPI_Comm				x_comm;
	static	int						is,ie;
	
	/*CREATE LOCAL PLANES*/
	static	MPI_Datatype			local_x_plane_red;
	static	MPI_Datatype			local_x_plane_blk;
	static	MPI_Datatype			x_plane;
	static	MPI_Datatype			y_plane;
	static	MPI_Datatype			z_plane;
	
	static	CMatrix2D				top_plane;
	static	CMatrix2D				bot_plane;
	static	CMatrix2D				front_plane;
	static	CMatrix2D				back_plane;
	static	CMatrix2D				left_plane;
	static	CMatrix2D				right_plane;
	
	/*SEND/RECV PARAMETERS*/
	static	int						tag;
	static	MPI_Status				status;
	static	MPI_Request				request;
};
/**********************************************************************************/

#endif // MPIINPUT_H
