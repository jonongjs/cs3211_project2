#include <mpi.h>  
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <math.h>
#define BUFSIZE 100000

int modexp(int g, int u,int p){
	int s = 1;
	//printf("\ng=%d p=%d u=%d",g,p,u);
	while (u!=0){
      if ((u & 1) !=0)
        s = (s * g ) % p;
      u = u >>  1;
      g = (g*g) % p;
    }
	//printf(" s=%d ",s);
	return s;
}

int main(int argc, char *argv[]){
	int factors[100000]={0};
	int gen[100000] = {0};
	int count[32]={0};
	int countgen[32] = {0};
	int numprocs,myid;
	int sqr,max_e_per_proc;
	int startval,i;
	int phi;
	int globalsum;
	const int root =0;
	int p;
	if (argc>=2){
		p = atoi(argv[1]);
	}
	else{
		printf("No value for p\n");
		return 0;
	}
	
	MPI_Status stat;
	MPI_Init(&argc,&argv); 
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	phi = p -1 ;
	//sqr = ceil(sqrt(phi));
	sqr = phi;
	while (sqr%numprocs) ++sqr;
	max_e_per_proc = sqr/numprocs;
	if (myid==root){
		printf("Process rank 0, p = %d, max_e_per_proc=%d\n", p, max_e_per_proc);
		startval = 2;
	}else{
		startval = max_e_per_proc*myid + 1;
		printf("Not root, rank %d startval= %d\n", myid, startval);
	}
	for (i=0;i<max_e_per_proc;++i){
		int curval = startval + i;
		if (curval >= phi) break;
		if ((phi % curval) == 0){
			factors[count[0]] = curval;
			//printf(" [%d] " , curval);
			++count[0];
		}
	}
	MPI_Allgather(factors, max_e_per_proc, MPI_INT, factors, max_e_per_proc, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(count, 1, MPI_INT, count , 1, MPI_INT, MPI_COMM_WORLD);
	
	i=0;
	int z=0;
	int flag = 0;
	//printf("The factors: ");
	for (z=0;z<max_e_per_proc;++z){
		int curval = startval + z;
		int upper = max_e_per_proc*(myid+1) +1;
		flag = 0;
		//printf("\ncurval = %d ,upper=%d ",curval,upper);
		if ((curval > p-1 )|| (curval >= upper)) break;
		for (i=0;i<numprocs;++i){
			int j=0;
			int stv = max_e_per_proc*i;
			for (j=0;j<count[i];j++){
				//printf("\ncurval=%d ", curval);
				if (modexp(curval, factors[stv+j],p)==1){
					flag = 1;
					break;
				}
			}
			if (flag) break;
		}
		if (!flag){
			//printf("\nfiltered curval=%d ", curval);
			gen[countgen[0]]=curval;
			++countgen[0];
		}
	}
	//}
	MPI_Reduce(&countgen[0], &globalsum, 1, MPI_INT, MPI_SUM , root, MPI_COMM_WORLD);
	MPI_Allgather(gen, max_e_per_proc, MPI_INT, gen, max_e_per_proc, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(countgen, 1, MPI_INT, countgen , 1, MPI_INT, MPI_COMM_WORLD);
	
	if (myid==root){
		printf("\ntotal number of generators:%d \n",globalsum);
		//printf("The generators: ");
		/*for (i=0;i<numprocs;++i){
			int j=0;
			int stv = max_e_per_proc*i;
			for (j=0;j<countgen[i];j++){
				printf("%d ",gen[stv+j]);
			}
		}*/
		
		printf("\n");
	}
	MPI_Finalize(); 
	return 0;
}
