/* FLPQ-HDF5 v1: independent per-quantity serial root writer.
 * MPI ranks evaluate bounded buffers; rank 0 writes HDF5 hyperslabs.
 * No dependency on h5py or PQViewer. 2026-10-20 F3 implementation.
 */
#include "flpq.h"
#include "Local_physical_quantities.h"
#include "FLPQ_HDF5.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/resource.h>

/* Limit peak temporary memory to a few hundred KiB per MPI rank. */
#define FLPQ_H5_MAX_BLOCK 8192
#define FLPQ_H5_MAX_COMPONENTS 12
#define FLPQ_H5_TAG_HEADER 7191
#define FLPQ_H5_TAG_VALUES 7192

static void h5_fail(const char *msg)
{
    int me;
    MPI_Comm_rank(mpi_comm_level1, &me);
    fprintf(stderr, "[FLPQ-HDF5 rank %d] %s\n", me, msg);
    fflush(stderr);
    MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
}
#define CHECK_H5(expr, msg) do { if ((expr) < 0) h5_fail(msg); } while (0)

static void astring(hid_t group, const char *key, const char *value)
{
    hid_t ds = H5Screate(H5S_SCALAR), type = H5Tcopy(H5T_C_S1), at;
    if (ds < 0 || type < 0) h5_fail("creating string attribute");
    CHECK_H5(H5Tset_size(type, strlen(value) + 1), "setting attr length");
    CHECK_H5(H5Tset_strpad(type, H5T_STR_NULLTERM), "setting attr padding");
    at = H5Acreate2(group, key, type, ds, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_H5(at, "creating string attr");
    CHECK_H5(H5Awrite(at, type, value), "writing string attr");
    CHECK_H5(H5Aclose(at), "closing string attr");
    CHECK_H5(H5Tclose(type), "closing string type");
    CHECK_H5(H5Sclose(ds), "closing attribute dataspace");
}

static void aint(hid_t group, const char *key, int value)
{
    hid_t space = H5Screate(H5S_SCALAR), at;
    CHECK_H5(space, "creating integer attribute space");
    at = H5Acreate2(group, key, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_H5(at, "creating integer attribute");
    CHECK_H5(H5Awrite(at, H5T_NATIVE_INT, &value), "writing integer attribute");
    CHECK_H5(H5Aclose(at), "closing integer attribute");
    CHECK_H5(H5Sclose(space), "closing integer dataspace");
}

static void dataset(hid_t group, const char *name, hid_t type,
                    hid_t native_type, int ndim, const hsize_t *dims, const void *ptr)
{
    hid_t space = H5Screate_simple(ndim, dims, NULL), d;
    CHECK_H5(space, "creating metadata space");
    d = H5Dcreate2(group, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_H5(d, "creating metadata dataset");
    if (ptr != NULL)
        CHECK_H5(H5Dwrite(d, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr),
                 "writing metadata dataset");
    CHECK_H5(H5Dclose(d), "closing metadata dataset");
    CHECK_H5(H5Sclose(space), "closing metadata space");
}

static hid_t field_create(hid_t grp, const char *name, int dims_count,
                          const hsize_t *dims, const hsize_t *chunks)
{
    hid_t sp = H5Screate_simple(dims_count, dims, NULL), pl, ds;
    CHECK_H5(sp, "creating HDF5 field dataspace");
    pl = H5Pcreate(H5P_DATASET_CREATE);
    CHECK_H5(pl, "creating HDF5 dataset property list");
    CHECK_H5(H5Pset_chunk(pl, dims_count, chunks), "setting chunk shape");
    if (FLPQ_HDF5_Compression == 1) {
        CHECK_H5(H5Pset_shuffle(pl), "setting shuffle filter");
        CHECK_H5(H5Pset_deflate(pl, FLPQ_HDF5_GzipLevel), "setting gzip compression");
    }
    ds = H5Dcreate2(grp, name, FLPQ_HDF5_Precision==32 ? H5T_IEEE_F32LE : H5T_IEEE_F64LE, sp,
                    H5P_DEFAULT, pl, H5P_DEFAULT);
    CHECK_H5(ds, "creating HDF5 field dataset");
    CHECK_H5(H5Pclose(pl), "closing HDF5 creation properties");
    CHECK_H5(H5Sclose(sp), "closing HDF5 field dataspace");
    return ds;
}

static void create_metadata(hid_t file, const char *quantity,
                            int rank, int components,
                            int nx, int ny, int nz)
{
    hid_t g, a, f;
    hsize_t s1[1], s2[2], dims4[4], chunks4[4], dims5[5], chunks5[5];
    long long dims[3];
    double origin[3], step[3][3], lattice[3][3];
    int periodic_known = 0, periodic[3] = {0,0,0};
    int *numbers = NULL;
    double *positions = NULL, *charges = NULL;
    int i,j,na = atomnum;
    const char *order = rank == 0 ? "scalar" :
                        rank == 1 ? "x y z" :
                        strcmp(quantity, "stress_diag") == 0 ?
                        "eigenvalues(3),eigenvectors(3x3)" :
                        "xx xy xz yx yy yz zx zy zz";

    astring(file, "format_id", "FLPQ-HDF5");
    aint(file, "schema_version", 1);
    astring(file, "producer", "FLPQ");
    astring(file, "quantity", quantity);
    astring(file, "coordinate_unit", "bohr");
    astring(file, "value_unit", "atomic_unit");
    aint(file, "field_rank", rank);
    astring(file, "precision", FLPQ_HDF5_Precision==32 ? "float32" : "float64");
    astring(file, "compression", FLPQ_HDF5_Compression==1 ? "gzip" : "none");
    aint(file, "gzip_level", FLPQ_HDF5_Compression==1 ? FLPQ_HDF5_GzipLevel : 0);
    astring(file, "component_basis", "cartesian");
    astring(file, "component_order", order);
    astring(file, "data_order", "component,n1,n2,n3");
    aint(file, "periodicity_known", periodic_known);
    aint(file, "source_spin_mode", SpinP_switch);

    dims[0]=nx;dims[1]=ny;dims[2]=nz;
    for(i=0;i<3;i++){
        origin[i]=Grid_Origin_N[0][i+1];
        for(j=0;j<3;j++){
            /* Ngrid1 is temporarily 1 in save_memory_3d; keep the full-grid step. */
            step[i][j]=(flag_save_memory_3d && i==0) ?
                tv[1][j+1]/(double)loop_Ngrid1 : gtv[i+1][j+1];
            lattice[i][j]=tv_ori[i+1][j+1];
        }
    }
    g=H5Gcreate2(file,"grid",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    CHECK_H5(g,"creating grid group");
    s1[0]=3;
    dataset(g,"shape",H5T_STD_I64LE,H5T_NATIVE_LLONG,1,s1,dims);
    dataset(g,"origin",H5T_IEEE_F64LE,H5T_NATIVE_DOUBLE,1,s1,origin);
    s2[0]=3;s2[1]=3;
    dataset(g,"step_vectors",H5T_IEEE_F64LE,H5T_NATIVE_DOUBLE,2,s2,step);
    dataset(g,"lattice_vectors",H5T_IEEE_F64LE,H5T_NATIVE_DOUBLE,2,s2,lattice);
    /* Periodicity not determined by FLPQ input; preserve unknown separately. */
    {
        unsigned char flags[3]={0,0,0};
        (void)periodic;
        dataset(g,"periodic",H5T_STD_U8LE,H5T_NATIVE_UCHAR,1,s1,flags);
    }
    CHECK_H5(H5Gclose(g),"closing grid group");

    if (na < 0) h5_fail("invalid atom count");
    numbers=(int*)malloc(sizeof(int)*(size_t)(na+1));
    positions=(double*)malloc(sizeof(double)*(size_t)(na+1)*3);
    charges=(double*)malloc(sizeof(double)*(size_t)(na+1));
    if (!numbers || !positions || !charges) h5_fail("atom metadata allocation failed");
    for (i=0;i<na;i++){
        int atom=i+1, spe=WhatSpecies[atom];
        numbers[i]=Spe_WhatAtom[spe];
        charges[i]=Spe_Core_Charge[spe]-InitN_USpin[atom]-InitN_DSpin[atom];
        for (j=0;j<3;j++) positions[3*i+j]=Gxyz[atom][j+1];
    }
    a=H5Gcreate2(file,"atoms",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    CHECK_H5(a,"creating atoms group");
    s1[0]=(hsize_t)na;
    s2[0]=(hsize_t)na;s2[1]=3;
    dataset(a,"atomic_numbers",H5T_STD_I32LE,H5T_NATIVE_INT,1,s1,numbers);
    dataset(a,"positions",H5T_IEEE_F64LE,H5T_NATIVE_DOUBLE,2,s2,positions);
    dataset(a,"cube_charges",H5T_IEEE_F64LE,H5T_NATIVE_DOUBLE,1,s1,charges);
    CHECK_H5(H5Gclose(a),"closing atoms group");
    free(numbers);free(positions);free(charges);

    f=H5Gcreate2(file,"field",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    CHECK_H5(f,"creating field group");
    chunks4[0]=1;chunks4[1]=1;
    chunks4[2]=(hsize_t)(ny<FLPQ_HDF5_ChunkY?ny:FLPQ_HDF5_ChunkY);
    chunks4[3]=(hsize_t)(nz<FLPQ_HDF5_ChunkZ?nz:FLPQ_HDF5_ChunkZ);
    if (strcmp(quantity,"stress_diag")==0){
        hid_t d;
        dims4[0]=3;dims4[1]=nx;dims4[2]=ny;dims4[3]=nz;
        d=field_create(f,"eigenvalues",4,dims4,chunks4);
        CHECK_H5(H5Dclose(d),"closing eigenvalues dataset");
        dims5[0]=3;dims5[1]=3;
        dims5[2]=nx;dims5[3]=ny;dims5[4]=nz;
        chunks5[0]=1;chunks5[1]=1;chunks5[2]=1;
        chunks5[3]=(hsize_t)(ny<FLPQ_HDF5_ChunkY?ny:FLPQ_HDF5_ChunkY);
        chunks5[4]=(hsize_t)(nz<FLPQ_HDF5_ChunkZ?nz:FLPQ_HDF5_ChunkZ);
        d=field_create(f,"eigenvectors",5,dims5,chunks5);
        CHECK_H5(H5Dclose(d),"closing eigenvectors dataset");
    }
    else{
        hid_t d;
        dims4[0]=components;dims4[1]=nx;dims4[2]=ny;dims4[3]=nz;
        d=field_create(f,"values",4,dims4,chunks4);
        CHECK_H5(H5Dclose(d),"closing field values dataset");
    }
    CHECK_H5(H5Gclose(f),"closing field group");
}

static void write_block(hid_t values, hid_t eigenvalues, hid_t eigenvectors,
                        const double *interleaved, int count, int components,
                        long long global_index, int nx, int ny, int nz)
{
    double *column = NULL;
    int allocated=FLPQ_HDF5_BlockPoints;
    int local=0, c,t;
    const long long per_plane=(long long)ny*nz, total=(long long)nx*ny*nz;
    column=(double*)malloc((size_t)allocated*sizeof(double));
    if (!column) h5_fail("allocating HDF5 output scratch buffer");
    if (global_index<0 || global_index+(long long)count>total)
        h5_fail("MPI block outside global grid");
    while(local<count){
        long long line_idx=global_index+(long long)local, rem;
        hsize_t start4[4], cnt4[4]={1,1,1,1};
        hsize_t start5[5], cnt5[5]={1,1,1,1,1};
        int run=nz-(int)(line_idx % nz);
        int rows=1, remaining=count-local;
        rem=line_idx%per_plane;
        /* Merge complete z-rows into a single contiguous HDF5 write.
         * Stop at a plane boundary; never cross n1 or mix components.
         * H5Sselect_hyperslab selects [1,1,rows,nz] or [1,1,1,rows,nz].
         */
        if(rem%nz==0 && remaining>=2*nz) {
            int available_rows=remaining/nz;
            int rows_left=ny-(int)(rem/nz);
            rows=available_rows<rows_left?available_rows:rows_left;
            if(rows>allocated/nz) rows=allocated/nz;
            if(rows>=2) run=rows*nz;
            else rows=1;
        }
        if(run>remaining)run=remaining;
        if(run>allocated) run=allocated; /* allow blocks smaller than Nz */
        for (c=0;c<components;c++){
            hid_t ds, file_space, mem_space;
            int dim;
            for(t=0;t<run;t++)column[t]=interleaved[(size_t)(local+t)*components+c];
            if (components==12){
                if (c<3){
                    ds=eigenvalues;dim=4;
                    start4[0]=c;start4[1]=line_idx/per_plane;
                    start4[2]=rem/nz;start4[3]=rem%nz;
                }
                else{
                    ds=eigenvectors;dim=5;
                    start5[0]=(c-3)/3;start5[1]=(c-3)%3;
                    start5[2]=line_idx/per_plane;
                    start5[3]=rem/nz;start5[4]=rem%nz;
                }
            }
            else{
                ds=values;dim=4;
                start4[0]=c;start4[1]=line_idx/per_plane;
                start4[2]=rem/nz;start4[3]=rem%nz;
            }
            file_space=H5Dget_space(ds);
            CHECK_H5(file_space,"getting field dataspace");
            if(dim==4){
                cnt4[2]=(hsize_t)rows;
                cnt4[3]=(hsize_t)(rows>1?nz:run);
                CHECK_H5(H5Sselect_hyperslab(file_space,H5S_SELECT_SET,start4,
                                             NULL,cnt4,NULL),"selecting field hyperslab");
            }
            else{
                cnt5[3]=(hsize_t)rows;
                cnt5[4]=(hsize_t)(rows>1?nz:run);
                CHECK_H5(H5Sselect_hyperslab(file_space,H5S_SELECT_SET,start5,
                                             NULL,cnt5,NULL),"selecting eigenvector hyperslab");
            }
            {
                hsize_t one[1]={(hsize_t)run};
                mem_space=H5Screate_simple(1,one,NULL);
            }
            CHECK_H5(mem_space,"creating HDF5 memory dataspace");
            CHECK_H5(H5Dwrite(ds,H5T_NATIVE_DOUBLE,mem_space,file_space,
                              H5P_DEFAULT,column),"writing field values");
            CHECK_H5(H5Sclose(mem_space),"closing memory dataspace");
            CHECK_H5(H5Sclose(file_space),"closing field dataspace");
        }
        local+=run;
    }
    free(column);
}

static void evaluate_range(int rank, int index, int start, int count,
                           int components, double *output)
{
    int i;
    for(i=0;i<count;i++){
        int c;
        FLPQ_LocalValue v;
        if (!flpq_evaluate_local_quantity(rank,index,start+i,&v)
             || v.components!=components)
            h5_fail("local physical quantity returned invalid components");
        for(c=0;c<components;c++){
            if (!isfinite(v.values[c])) h5_fail("nonfinite local quantity");
            output[(size_t)i*components+c]=v.values[c];
        }
    }
}

void FLPQ_HDF5_write(const char *name, int rank, int index)
{
    int me, nprocs, p, nx, ny, nz, components, block;
    double wall_start = MPI_Wtime(), local_eval = 0.0, max_eval = 0.0;
    double io_seconds = 0.0, max_wall = 0.0, local_wall;
    double sample_start;
    long max_rss = 0, max_rss_all = 0;
    struct rusage usage;
    long long measured_bytes = -1;
    long long expected=0, all=(long long)Ngrid1*Ngrid2*Ngrid3;
    long long start, count, meta[2];
    double *buffer=NULL;
    char partial[YOUSO10], final[YOUSO10];
    const char *quantity=(name[0]=='.'?name+1:name);
    hid_t file=-1, values=-1, eig=-1, vec=-1;
    MPI_Status status;
    MPI_Comm_rank(mpi_comm_level1,&me);
    MPI_Comm_size(mpi_comm_level1,&nprocs);

    nx=(flag_save_memory_3d ? loop_Ngrid1:Ngrid1);
    ny=Ngrid2;nz=Ngrid3;
    components=(rank==0 ? 1:rank==1?3:
                (rank==2 && index==3)?12:9);
    if(nx<=0 || ny<=0 || nz<=0 || components>12 ||
       all<=0 || all>(long long)INT_MAX || My_NumGridB_AB<0)
        h5_fail("invalid grid, count or field components");

    /* FLPQ's distributed B-grid is a contiguous set of z-fastest rows. */
    start=((long long)me*(Ngrid1*Ngrid2)+nprocs-1)/nprocs*(long long)Ngrid3;
    count=My_NumGridB_AB;
    block=FLPQ_HDF5_BlockPoints;
    buffer=(double*)malloc((size_t)block*components*sizeof(double));
    if (!buffer)h5_fail("allocating bounded HDF5 send buffer");

    if (snprintf(partial,sizeof(partial),"%s/LPQ.%s.partial.h5",
                 DirnameLPQ,quantity)>=(int)sizeof(partial) ||
        snprintf(final,sizeof(final),"%s/LPQ.%s.h5",
                 DirnameLPQ,quantity)>=(int)sizeof(final))
        h5_fail("HDF5 file name exceeds buffer");

    if(me==Host_ID){
        sample_start=MPI_Wtime();
        if(N_loop_Ngrid1==0){
            struct stat st;
            if(lstat(partial,&st)==0)
                h5_fail("existing partial HDF5 found; remove/resolve it manually");
            if(lstat(final,&st)==0) {
                if(FLPQ_HDF5_OnExisting==0)
                    h5_fail("existing HDF5 output found; set Output.HDF5.OnExisting replace explicitly");
                if(!S_ISREG(st.st_mode))
                    h5_fail("existing HDF5 target is not a regular file");
            }
            file=H5Fcreate(partial,H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
            CHECK_H5(file,"creating new partial HDF5 file");
            create_metadata(file,quantity,rank,components,nx,ny,nz);
        }
        else{
            file=H5Fopen(partial,H5F_ACC_RDWR,H5P_DEFAULT);
            CHECK_H5(file,"opening existing partial HDF5 for the next slice");
        }
        if(components==12){
            eig=H5Dopen2(file,"/field/eigenvalues",H5P_DEFAULT);
            vec=H5Dopen2(file,"/field/eigenvectors",H5P_DEFAULT);
            CHECK_H5(eig,"opening eigenvalue dataset");
            CHECK_H5(vec,"opening eigenvector dataset");
        }
        else{
            values=H5Dopen2(file,"/field/values",H5P_DEFAULT);
            CHECK_H5(values,"opening quantity dataset");
        }
        io_seconds += MPI_Wtime()-sample_start;
    }

    /* Point ranges are exchanged without a global-sized gather on rank 0. */
    if(me!=Host_ID){
        meta[0]=start;meta[1]=count;
        MPI_Send(meta,2,MPI_LONG_LONG_INT,Host_ID,FLPQ_H5_TAG_HEADER,mpi_comm_level1);
        for(p=0;p<count;p+=FLPQ_HDF5_BlockPoints){
            int n=(count-p<FLPQ_HDF5_BlockPoints?count-p:FLPQ_HDF5_BlockPoints);
            sample_start=MPI_Wtime();
            evaluate_range(rank,index,p,n,components,buffer);
            local_eval+=MPI_Wtime()-sample_start;
            MPI_Send(buffer,n*components,MPI_DOUBLE,Host_ID,
                     FLPQ_H5_TAG_VALUES,mpi_comm_level1);
        }
    }
    else{
        for(p=0;p<nprocs;p++){
            long long begin, nlocal, j;
            if(p==Host_ID){begin=start;nlocal=count;}
            else{
                MPI_Recv(meta,2,MPI_LONG_LONG_INT,p,FLPQ_H5_TAG_HEADER,
                         mpi_comm_level1,&status);
                begin=meta[0];nlocal=meta[1];
            }
            if(begin!=expected || nlocal<0 || begin+nlocal>all)
                h5_fail("MPI grid rank partition is non-contiguous or incomplete");
            for(j=0;j<nlocal;j+=FLPQ_HDF5_BlockPoints){
                int n=(nlocal-j<FLPQ_HDF5_BlockPoints?(int)(nlocal-j):FLPQ_HDF5_BlockPoints);
                if(p==Host_ID) {
                    sample_start=MPI_Wtime();
                    evaluate_range(rank,index,(int)j,n,components,buffer);
                    local_eval+=MPI_Wtime()-sample_start;
                }
                else MPI_Recv(buffer,n*components,MPI_DOUBLE,p,
                              FLPQ_H5_TAG_VALUES,mpi_comm_level1,&status);
                sample_start=MPI_Wtime();
                write_block(values,eig,vec,buffer,n,components,
                            (long long)N_loop_Ngrid1*all+begin+j,
                            nx,ny,nz);
                io_seconds+=MPI_Wtime()-sample_start;
            }
            expected=begin+nlocal;
        }
        if(expected!=all)h5_fail("MPI grid rank partition does not cover the full slice");
        sample_start=MPI_Wtime();
        if(values>=0)CHECK_H5(H5Dclose(values),"closing field dataset");
        if(eig>=0)CHECK_H5(H5Dclose(eig),"closing eigenvalue dataset");
        if(vec>=0)CHECK_H5(H5Dclose(vec),"closing eigenvector dataset");
        CHECK_H5(H5Fclose(file),"closing partial HDF5 file");
        if(N_loop_Ngrid1==loop_Ngrid1-1){
            if(rename(partial,final)!=0)
                h5_fail("renaming completed HDF5 output");
        }
        io_seconds+=MPI_Wtime()-sample_start;
        if(N_loop_Ngrid1==loop_Ngrid1-1){
            struct stat st;
            if(stat(final,&st)==0) measured_bytes=(long long)st.st_size;
        }
    }
    free(buffer);
    /* Report a per-quantity full-slice timing; HDF5 I/O alone excludes
       numerical evaluation and MPI waiting. ru_maxrss is process high-water. */
    if(FLPQ_HDF5_LogPerformance) {
        getrusage(RUSAGE_SELF,&usage);
        max_rss=usage.ru_maxrss;
        local_wall=MPI_Wtime()-wall_start;
        MPI_Reduce(&local_eval,&max_eval,1,MPI_DOUBLE,MPI_MAX,Host_ID,mpi_comm_level1);
        MPI_Reduce(&local_wall,&max_wall,1,MPI_DOUBLE,MPI_MAX,Host_ID,mpi_comm_level1);
        MPI_Reduce(&max_rss,&max_rss_all,1,MPI_LONG,MPI_MAX,Host_ID,mpi_comm_level1);
        if(me==Host_ID) {
            printf("FLPQ_HDF5_PERF quantity=%s slice=%d/%d ranks=%d nx=%d ny=%d nz=%d "
                   "block_points=%d chunk_y=%d chunk_z=%d precision=%d gzip=%d "
                   "gzip_level=%d io_seconds=%.6f wall_seconds=%.6f "
                   "eval_max_seconds=%.6f max_rss_kib=%ld file_bytes=%lld\n",
                   quantity,N_loop_Ngrid1+1,loop_Ngrid1,nprocs,nx,ny,nz,
                   block,FLPQ_HDF5_ChunkY,FLPQ_HDF5_ChunkZ,
                   FLPQ_HDF5_Precision,FLPQ_HDF5_Compression,
                   FLPQ_HDF5_GzipLevel,io_seconds,max_wall,max_eval,
                   max_rss_all,measured_bytes);
            fflush(stdout);
        }
    }
    MPI_Barrier(mpi_comm_level1);
}
