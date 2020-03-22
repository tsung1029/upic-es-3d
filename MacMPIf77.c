/* This is a Fortran callable interface library to the C version of a
   Partial MPI library based on the Program-to-Program Communications
   ToolBox in the Macintosh OS.  No local buffering of messages is
   implemented, so that all messages must be received in the order sent,
   and receives with wildcard sources are not supported.
   the following subroutines are implemented:
   MPI_Init, MPI_Finalize, MPI_Send, MPI_Recv, MPI_Isend, MPI_Irecv
   MPI_Test, MPI_Wait, MPI_Sendrecv, MPI_Ssend, MPI_Issend, MPI_Waitall
   MPI_Waitany, MPI_Get_count, MPI_Initialized, MPI_Comm_size
   MPI_Comm_rank, MPI_Comm_dup, MPI_Comm_split, MPI_Comm_free
   MPI_Cart_create, MPI_Cart_coords, MPI_Cart_get, MPI_Cart_shift
   MPI_Cart_rank, MPI_Cart_sub, MPI_Dims_create
   MPI_Bcast, MPI_Barrier, MPI_Reduce, MPI_Scan
   MPI_Allreduce, MPI_Gather, MPI_Allgather, MPI_Scatter, MPI_Alltoall
   MPI_Gatherv, MPI_Allgatherv, MPI_Scatterv, MPI_Alltoallv
   MPI_Reduce_scatter, MPI_Abort, MPI_Wtime, MPI_Wtick, MPI_Type_extent
   MPI_Request_free, MPI_Get_processor_name, MPI_Errhandler_set
   written by viktor k. decyk, ucla
   copyright 1998, regents of the university of california
   update: april 1, 2004                                                 */

#include <stdlib.h>
#include "mpi.h"

#define MPI_REAL                19
#define MPI_DOUBLE_PRECISION    20
#define MPI_COMPLEX             22
#define MPI_DOUBLE_COMPLEX      23

void mpi_init_(int *ierror) {
   *ierror = MPI_Init(NULL,NULL);
   return;
}

void mpi_finalize_(int *ierror) {
   *ierror = MPI_Finalize();
   return;
}

void mpi_send_(void* buf, int *count, MPI_Datatype *datatype, int *dest,
             int *tag, MPI_Comm *comm, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Send(buf,lcount,ldatatype,*dest,*tag,*comm);
   return;
}

void mpi_recv_(void* buf, int *count, MPI_Datatype *datatype, int *source,
             int *tag, MPI_Comm *comm, MPI_Status *status, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Recv(buf,lcount,ldatatype,*source,*tag,*comm,status);
   return;
}

void mpi_isend_(void* buf, int *count, MPI_Datatype *datatype, int *dest,
              int *tag, MPI_Comm *comm, MPI_Request *request, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Isend(buf,lcount,ldatatype,*dest,*tag,*comm,request);
   return;
}

void mpi_irecv_(void* buf, int *count, MPI_Datatype *datatype, int *source,
              int *tag, MPI_Comm *comm, MPI_Request *request, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Irecv(buf,lcount,ldatatype,*source,*tag,*comm,request);
   return;
}

void mpi_test_(MPI_Request *request, int *flag, MPI_Status *status,
             int *ierror) {
   *ierror = MPI_Test(request,flag,status);
   return;
}

void mpi_wait_(MPI_Request *request, MPI_Status *status, int *ierror) {
   *ierror = MPI_Wait(request,status);
   return;
}

void mpi_request_free_(MPI_Request *request, int *ierror) {
   *ierror = MPI_Request_free(request);
   return;
}

void mpi_sendrecv_(void* sendbuf, int *sendcount, MPI_Datatype *sendtype,
                 int *dest, int *sendtag, void* recvbuf, int *recvcount,
                 MPI_Datatype *recvtype, int *source, int *recvtag,
                 MPI_Comm *comm, MPI_Status *status, int *ierror) {
   int lsendcount, lrecvcount;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   *ierror = MPI_Sendrecv(sendbuf,lsendcount,lsendtype,*dest,*sendtag,
                          recvbuf,lrecvcount,lrecvtype,*source,*recvtag,
                          *comm,status);
   return;
}

void mpi_ssend_(void* buf, int *count, MPI_Datatype *datatype, int *dest,
              int *tag, MPI_Comm *comm, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Ssend(buf,lcount,ldatatype,*dest,*tag,*comm);
   return;
}

void mpi_issend_(void* buf, int *count, MPI_Datatype *datatype, int *dest,
               int *tag, MPI_Comm *comm, MPI_Request *request, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Issend(buf,lcount,ldatatype,*dest,*tag,*comm,request);
   return;
}

void mpi_waitall_(int *count, MPI_Request *array_of_requests,
                MPI_Status *array_of_statuses, int *ierror) {
   *ierror = MPI_Waitall(*count,array_of_requests,array_of_statuses);
   return;
}

void mpi_waitany_(int *count, MPI_Request *array_of_requests,
                int *index, MPI_Status *status, int *ierror) {
   *ierror = MPI_Waitany(*count,array_of_requests,index,status);
   return;
}

void mpi_get_count_(MPI_Status *status, MPI_Datatype *datatype,
                  int *count, int *ierror) {
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX)
      ldatatype = MPI_REAL;
   else if (*datatype==MPI_DOUBLE_COMPLEX)
      ldatatype = MPI_DOUBLE_PRECISION;
   else
      ldatatype = *datatype;
   *ierror = MPI_Get_count(status,ldatatype,count);
   if (*datatype==MPI_COMPLEX)
      *count = *count/2;
   else if (*datatype==MPI_DOUBLE_COMPLEX)
      *count = *count/2;
   return;
}

void mpi_initialized_(int *flag, int *ierror) {
   *ierror = MPI_Initialized(flag);
   return;
}

void mpi_comm_size_(MPI_Comm *comm, int *size, int *ierror) {
   *ierror = MPI_Comm_size(*comm,size);
   return;
}

void mpi_comm_rank_(MPI_Comm *comm, int *rank, int *ierror) {
   *ierror = MPI_Comm_rank(*comm,rank);
   return;
}

void mpi_comm_dup_(MPI_Comm *comm, MPI_Comm *newcomm, int *ierror) {
   *ierror = MPI_Comm_dup(*comm,newcomm);
   return;
}

void mpi_comm_split_(MPI_Comm *comm, int *color, int *key, MPI_Comm *newcomm,
                    int *ierror) {
   *ierror = MPI_Comm_split(*comm,*color,*key,newcomm);
   return;
}

void mpi_comm_free_(MPI_Comm *comm, int *ierror) {
   *ierror = MPI_Comm_free(comm);
   return;
}

void mpi_cart_create_(MPI_Comm *comm_old, int *ndims, int *dims, int *periods,
                     int *reorder, MPI_Comm *comm_cart,int *ierror) {
   *ierror = MPI_Cart_create(*comm_old,*ndims,dims,periods,*reorder,comm_cart);
   return;
}

void mpi_cart_coords_(MPI_Comm *comm, int *rank, int *maxdims, int *coords,
                     int *ierror) {
   *ierror = MPI_Cart_coords(*comm,*rank,*maxdims,coords);
   return;
}

void mpi_cart_get_(MPI_Comm *comm, int *maxdims, int *dims, int *periods,
                 int *coords, int *ierror) {
   *ierror = MPI_Cart_get(*comm,*maxdims,dims,periods,coords);
   return;
}

void mpi_cart_shift_(MPI_Comm *comm, int *direction, int *disp, int *rank_source,
                    int *rank_dest, int *ierror) {
   *ierror = MPI_Cart_shift(*comm,*direction,*disp,rank_source,rank_dest);
   return;
}

void mpi_cart_rank_(MPI_Comm *comm, int *coords, int *rank, int *ierror) {
   *ierror = MPI_Cart_rank(*comm,coords,rank);
   return;
}

void mpi_cart_sub_(MPI_Comm *comm, int *remain_dims, int *newcomm, int *ierror) {
   *ierror = MPI_Cart_sub(*comm,remain_dims,newcomm);
   return;
}

void mpi_dims_create_(int *nnodes, int *ndims, int *dims, int *ierror) {
   *ierror = MPI_Dims_create(*nnodes,*ndims,dims);
   return;
}

void mpi_bcast_(void* buffer, int *count, MPI_Datatype *datatype,
              int *root, MPI_Comm *comm, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_REAL;
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      lcount = 2*(*count);
      ldatatype = MPI_DOUBLE_PRECISION;
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Bcast(buffer,lcount,ldatatype,*root,*comm);
   return;
}

void mpi_barrier_(MPI_Comm *comm, int *ierror) {
   *ierror = MPI_Barrier(*comm);
   return;
}

void mpi_reduce_(void* sendbuf, void* recvbuf, int *count,
               MPI_Datatype *datatype, MPI_Op *op, int *root,
               MPI_Comm *comm, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_REAL;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_DOUBLE_PRECISION;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Reduce(sendbuf,recvbuf,lcount,ldatatype,*op,*root,*comm);
   return;
}

void mpi_scan_(void* sendbuf, void* recvbuf, int *count,
             MPI_Datatype *datatype, MPI_Op *op, MPI_Comm *comm, int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_REAL;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_DOUBLE_PRECISION;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Scan(sendbuf,recvbuf,lcount,ldatatype,*op,*comm);
   return;
}

void mpi_allreduce_(void* sendbuf, void* recvbuf, int *count,
                  MPI_Datatype *datatype, MPI_Op *op, MPI_Comm *comm,
                  int *ierror) {
   int lcount;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_REAL;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      if (*op==MPI_SUM) {
         lcount = 2*(*count);
         ldatatype = MPI_DOUBLE_PRECISION;
      }
      else {
         *ierror = 7;
         return;
      }
   }
   else {
      lcount = *count;
      ldatatype = *datatype;
   }
   *ierror = MPI_Allreduce(sendbuf,recvbuf,lcount,ldatatype,*op,*comm);
   return;
}

void mpi_gather_(void* sendbuf, int *sendcount, MPI_Datatype *sendtype,
               void* recvbuf, int *recvcount, MPI_Datatype *recvtype,
               int *root, MPI_Comm *comm, int *ierror) {
   int lsendcount, lrecvcount;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   *ierror = MPI_Gather(sendbuf,lsendcount,lsendtype,recvbuf,lrecvcount,
                       lrecvtype,*root,*comm);
   return;
}

void mpi_allgather_(void* sendbuf, int *sendcount,
                  MPI_Datatype *sendtype, void* recvbuf, int *recvcount,
                  MPI_Datatype *recvtype, MPI_Comm *comm, int *ierror) {
   int lsendcount, lrecvcount;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   *ierror = MPI_Allgather(sendbuf,lsendcount,lsendtype,recvbuf,lrecvcount,
                           lrecvtype,*comm);
   return;
}

void mpi_scatter_(void* sendbuf, int *sendcount, MPI_Datatype *sendtype,
                void* recvbuf, int *recvcount, MPI_Datatype *recvtype,
                int *root, MPI_Comm *comm, int *ierror) {
   int lsendcount, lrecvcount;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   *ierror = MPI_Scatter(sendbuf,lsendcount,lsendtype,recvbuf,lrecvcount,
                         lrecvtype,*root,*comm);
   return;
}

void mpi_alltoall_(void* sendbuf, int *sendcount, MPI_Datatype *sendtype,
                 void* recvbuf, int *recvcount, MPI_Datatype *recvtype,
                 MPI_Comm *comm, int *ierror) {
   int lsendcount, lrecvcount;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   *ierror = MPI_Alltoall(sendbuf,lsendcount,lsendtype,recvbuf,lrecvcount,
                          lrecvtype,*comm);
      return;
}

void mpi_gatherv_(void* sendbuf, int *sendcount, MPI_Datatype *sendtype,
                void* recvbuf, int *recvcounts, int *displs,
                MPI_Datatype *recvtype, int *root, MPI_Comm *comm,
                int *ierror) {
   int lsendcount, nproc, i;
   int *lrecvcounts, *ldispls;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX)
      lrecvtype = MPI_REAL;
   else if (*recvtype==MPI_DOUBLE_COMPLEX)
      lrecvtype = MPI_DOUBLE_PRECISION;
   else
      lrecvtype = *recvtype;
   if (lrecvtype != *recvtype) {
      i = MPI_Comm_size(*comm,&nproc);
      lrecvcounts = malloc(nproc);
      ldispls = malloc(nproc);
      if ((lrecvcounts==NULL) || (ldispls==NULL)) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         lrecvcounts[i] = 2*recvcounts[i];
         ldispls[i] = 2*displs[i];
      }
   }
   else {
      lrecvcounts = recvcounts;
      ldispls = displs;
   }
   *ierror = MPI_Gatherv(sendbuf,lsendcount,lsendtype,recvbuf,lrecvcounts,
                         ldispls,lrecvtype,*root,*comm);
   if (lrecvtype != *recvtype) {
      free(lrecvcounts);
      free(ldispls);
   }
   return;
}

void mpi_allgatherv_(void* sendbuf, int *sendcount,
                   MPI_Datatype *sendtype, void* recvbuf, int *recvcounts,
                   int *displs, MPI_Datatype *recvtype, MPI_Comm *comm,
                   int *ierror) {
   int lsendcount, nproc, i;
   int *lrecvcounts, *ldispls;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_REAL;
   }
   else if (*sendtype==MPI_DOUBLE_COMPLEX) {
      lsendcount = 2*(*sendcount);
      lsendtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lsendcount = *sendcount;
      lsendtype = *sendtype;
   }
   if (*recvtype==MPI_COMPLEX)
      lrecvtype = MPI_REAL;
   else if (*recvtype==MPI_DOUBLE_COMPLEX)
      lrecvtype = MPI_DOUBLE_PRECISION;
   else
      lrecvtype = *recvtype;
   if (lrecvtype != *recvtype) {
      i = MPI_Comm_size(*comm,&nproc);
      lrecvcounts = malloc(nproc);
      ldispls = malloc(nproc);
      if ((lrecvcounts==NULL) || (ldispls==NULL)) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         lrecvcounts[i] = 2*recvcounts[i];
         ldispls[i] = 2*displs[i];
      }
   }
   else {
      lrecvcounts = recvcounts;
      ldispls = displs;
   }
   *ierror = MPI_Allgatherv(sendbuf,lsendcount,lsendtype,recvbuf,
                            lrecvcounts,ldispls,lrecvtype,*comm);
   if (lrecvtype != *recvtype) {
      free(lrecvcounts);
      free(ldispls);
   }
   return;
}

void mpi_scatterv_(void* sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype *sendtype, void* recvbuf, int *recvcount,
                 MPI_Datatype *recvtype, int *root, MPI_Comm *comm,
                 int *ierror) {
   int nproc, lrecvcount, i;
   int *lsendcounts, *ldispls;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX)
      lsendtype = MPI_REAL;
   else if (*sendtype==MPI_DOUBLE_COMPLEX)
      lsendtype = MPI_DOUBLE_PRECISION;
   else
      lsendtype = *sendtype;
   if (*recvtype==MPI_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_REAL;
   }
   else if (*recvtype==MPI_DOUBLE_COMPLEX) {
      lrecvcount = 2*(*recvcount);
      lrecvtype = MPI_DOUBLE_PRECISION;
   }
   else {
      lrecvcount = *recvcount;
      lrecvtype = *recvtype;
   }
   if (lsendtype != *sendtype) {
      i = MPI_Comm_size(*comm,&nproc);
      lsendcounts = malloc(nproc);
      ldispls = malloc(nproc);
      if ((lsendcounts==NULL) || (ldispls==NULL)) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         lsendcounts[i] = 2*sendcounts[i];
         ldispls[i] = 2*displs[i];
      }
   }
   else {
      lsendcounts = sendcounts;
      ldispls = displs;
   }
   *ierror = MPI_Scatterv(sendbuf,lsendcounts,ldispls,lsendtype,recvbuf,
                         lrecvcount,lrecvtype,*root,*comm);
   if (lsendtype != *sendtype) {
      free(lsendcounts);
      free(ldispls);
   }
   return;
}

void mpi_alltoallv_(void* sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype *sendtype, void* recvbuf, int *recvcounts,
                  int *rdispls, MPI_Datatype *recvtype, MPI_Comm *comm,
                  int *ierror) {
   int nproc, i;
   int *lsendcounts, *lsdispls, *lrecvcounts, *lrdispls;
   MPI_Datatype lsendtype, lrecvtype;
   if (*sendtype==MPI_COMPLEX)
      lsendtype = MPI_REAL;
   else if (*sendtype==MPI_DOUBLE_COMPLEX)
      lsendtype = MPI_DOUBLE_PRECISION;
   else
      lsendtype = *sendtype;
   if (*recvtype==MPI_COMPLEX)
      lrecvtype = MPI_REAL;
   else if (*recvtype==MPI_DOUBLE_COMPLEX)
      lrecvtype = MPI_DOUBLE_PRECISION;
   else
      lrecvtype = *recvtype;
   if (lsendtype != *sendtype) {
      i = MPI_Comm_size(*comm,&nproc);
      lsendcounts = malloc(nproc);
      lsdispls = malloc(nproc);
      if ((lsendcounts==NULL) || (lsdispls==NULL)) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         sendcounts[i] = 2*sendcounts[i];
         lsdispls[i] = 2*sdispls[i];
      }
   }
   else {
      lsendcounts = sendcounts;
      lsdispls = sdispls;
   }
   if (lrecvtype != *recvtype) {
      i = MPI_Comm_size(*comm,&nproc);
      lrecvcounts = malloc(nproc);
      lrdispls = malloc(nproc);
      if ((lrecvcounts==NULL) || (lrdispls==NULL)) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         lrecvcounts[i] = 2*recvcounts[i];
         lrdispls[i] = 2*rdispls[i];
      }
   }
   else {
      lrecvcounts = recvcounts;
      lrdispls = rdispls;
   }
   *ierror = MPI_Alltoallv(sendbuf,lsendcounts,lsdispls,lsendtype,recvbuf,
                           lrecvcounts,lrdispls,lrecvtype,*comm);
   if (lsendtype != *sendtype) {
      free(lsendcounts);
      free(lsdispls);
   }
   if (lrecvtype != *recvtype) {
      free(recvcounts);
      free(lrdispls);
   }
   return;
}

void mpi_reduce_scatter_(void* sendbuf, void* recvbuf, int *recvcounts,
                       MPI_Datatype *datatype, MPI_Op *op, MPI_Comm *comm,
                       int *ierror) {
   int nproc, i;
   int *lrecvcounts;
   MPI_Datatype ldatatype;
   if (*datatype==MPI_COMPLEX) {
      if (*op==MPI_SUM)
         ldatatype = MPI_REAL;
      else {
         *ierror = 7;
         return;
      }
   }
   else if (*datatype==MPI_DOUBLE_COMPLEX) {
      if (*op==MPI_SUM)
         ldatatype = MPI_DOUBLE_PRECISION;
      else {
         *ierror = 7;
         return;
      }
   }
   else
      ldatatype = *datatype;
   if (ldatatype != *datatype) {
      i = MPI_Comm_size(*comm,&nproc);
      lrecvcounts = malloc(nproc);
      if (lrecvcounts==NULL) {
         *ierror = 21;
         return;
      }
      for (i = 0; i < nproc; i++) {
         lrecvcounts[i] = 2*recvcounts[i];
      }
   }
   else
      lrecvcounts = recvcounts;
   *ierror = MPI_Reduce_scatter(sendbuf,recvbuf,lrecvcounts,ldatatype,
                                *op,*comm);
   if (ldatatype != *datatype)
      free(lrecvcounts);
   return;
}

void mpi_abort_(MPI_Comm *comm, int *errorcode, int *ierror) {
   *ierror = MPI_Abort(*comm,*errorcode);
   return;
}

double mpi_wtime_(void) {
   return MPI_Wtime();
}

double mpi_wtick_(void) {
   return MPI_Wtick();
}

void mpi_type_extent_(MPI_Datatype *datatype, MPI_Aint *extent, int *ierror) {
   *ierror = MPI_Type_extent(*datatype,extent);
   return;
}

void mpi_errhandler_set_(MPI_Comm *comm, MPI_Errhandler *errhandler,
                         int *ierror) {
   *ierror = MPI_Errhandler_set(*comm,*errhandler);
   return;
}

void mpi_get_processor_name_(char *name, int *resultlen, int *ierror) {
   *ierror = MPI_Get_processor_name(name,resultlen);
   return;
}

void logname_(char* name, int n) {
   int i;
   char lname[37];
   if (n < 0)
      n = 0;
   else if (n > 36)
      n = 36;
   for (i = 0; i < n; i++) {
      lname[i] = name[i];
   }
   lname[n] = '\0';
   Logname(lname);
   return;
}

void set_mon_(int *monval) {
   Set_Mon(*monval);
   return;
}

int get_mon_() {
   return Get_Mon();
}
