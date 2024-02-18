#include <stdio.h>
#include <time.h>
#include <mpi.h>
#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255

struct complex{
  double real;
  double imag;
};
int cal_pixel(struct complex c) {
    double z_real = 0;
    double z_imag = 0;

    double z_real2, z_imag2, lengthsq;

    int iter = 0;
    do {
        z_real2 = z_real * z_real;
        z_imag2 = z_imag * z_imag;

        z_imag = 2 * z_real * z_imag + c.imag;
        z_real = z_real2 - z_imag2 + c.real;
        lengthsq =  z_real2 + z_imag2;
        iter++;
    }
    while ((iter < MAX_ITER) && (lengthsq < 4.0));

        return iter;

}
void save_pgm(const char *filename, int image[HEIGHT][WIDTH]) {
    FILE* pgmimg; 
    int temp;
    pgmimg = fopen(filename, "wb"); 
    fprintf(pgmimg, "P2\n"); // Writing Magic Number to the File   
    fprintf(pgmimg, "%d %d\n", WIDTH, HEIGHT);  // Writing Width and Height
    fprintf(pgmimg, "255\n");  // Writing the maximum gray value 
    int count = 0; 
    
    for (int i = 0; i < HEIGHT; i++) { 
        for (int j = 0; j < WIDTH; j++) { 
            temp = image[i][j]; 
            fprintf(pgmimg, "%d ", temp); // Writing the gray values in the 2D array to the file 
        } 
        fprintf(pgmimg, "\n"); 
    } 
    fclose(pgmimg); 
} 

int main(int argc, char const *argv[])
{
	int image[HEIGHT][WIDTH];
    struct complex c;
    MPI_Init(NULL,NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int i,j,l,row,s=HEIGHT/(world_size-1),slave;
    double time;
    
    printf("hello , i am %d \n", world_rank);
    if(world_rank==0){
        int colors[s][WIDTH+1];
        clock_t start_time = clock(); 
        for(i=1,row=0;i<world_size;i++,row+=s){
            MPI_Send(&row,1,MPI_INT,i,0, MPI_COMM_WORLD);
        }
        for(i=1;i<world_size;i++){
            MPI_Recv(&colors,s*(WIDTH+1)*sizeof(int),MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            slave=colors[0][WIDTH];
            for(j=slave, l=0;j<s+slave;j++,l++){
                //printf("1");
                for(int k=0;k<WIDTH;k++){
                    //printf("2");
                    image[j][k]=colors[l][k];
                }
            }
        }
        printf("\n3\n");
        clock_t end_time = clock();
        save_pgm("mandelbrot_dynamic.pgm", image);
        time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
        printf("The average execution time is: %f ms", time*1000);
    }
    else{
        int colors[s][WIDTH+1];
        MPI_Recv(&row,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for (i = row,l=0; i < row+s; i++,l++) {
            for (j = 0; j < WIDTH; j++) {
                c.real = (j - WIDTH / 2.0) * 4.0 / WIDTH;
                c.imag = (i - HEIGHT / 2.0) * 4.0 / HEIGHT;
                colors[l][j] = cal_pixel(c);
            }
        }
        printf("\ndone%d\n", world_rank);
        colors[0][WIDTH]=row;
        MPI_Send(&colors,s*(WIDTH+1)*sizeof(int),MPI_INT,0,0,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();

	return 0;
}