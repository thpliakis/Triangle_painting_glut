#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

//float canvas[512][512];

#define ALPHA 0.5
typedef struct 
{
	GLint v0[2];
	GLint v1[2];
	GLint v2[2];
}triangle;

typedef struct 
{
	GLfloat c0[3];
	GLfloat c1[3];
	GLfloat c2[3];
}color;

typedef struct
{
	int verts2d[4999][2];
	float vcolors[4999][3];
	int faces[10000][3];
	float depth[4999][1];
}data;

data d;

typedef struct 
{
	int x;
	int y;
	float rgb[3];
}point;


float interpolate_color(point p1, point p2, point p){
	float r1  = sqrt(pow(p2.x - p.x,2) + pow(p2.y - p.y,2));
	float r2  = sqrt(pow(p.x - p1.x,2) + pow(p.y - p1.y,2));
	float percent = r1 / (r1+r2);
	for(int i=0; i<3; i++){
		p.rgb[i] = percent * p1.rgb[i] + (1 - percent) * p2.rgb[i];
	}
}


void read_data_float( int M, int N, float m[M][N],  FILE *file){

	char line[250];

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			if(fgets(line, sizeof(line), file)){
				m[i][j] = atof(line);
			}else {
				break;
			}
		}
	}
}

void read_data_int( int M, int N, int m[M][N],  FILE *file){

	char line[250];

	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			if(fgets(line, sizeof(line), file)){
				m[i][j] = (int)atof(line);
			}else {
				break;
			}
		}
	}
}

void swap_int(int *a, int k, int l) {
	int temp = a[k];
	a[k] = a[l];
	a[l] = temp;
}

void swap_float(float *a, int k, int l) {
	float temp = a[k];
	a[k] = a[l];
	a[l] = temp;
}


/* partition -- in-place update of elements */
int partition2(float *a,int *o, int n) {
	float pivot = a[n-1];
	int i = 0;
	
	for (int j = 0; j < n - 1; j++){
		if (a[j] >= pivot){ 
		swap_int(o,i,j);
		swap_float(a,i++,j);
		}
	}

	swap_float(a, i, n - 1);
	swap_int(o, i, n - 1);
	return (i);
}

/* qsortseq -- Entry point for QuickSort */
void q_sort(float *a,int *o, int n) {
		
	if (n > 1) {
		//printf("n = %d \n",n);
		// if(l=3)
		// 	printf("a[%d] = %d\n",l,o[l]);
		// 	printf("a[%d] = %d\n",r,o[r]);
		int p = partition2(a,o, n);
		q_sort(a,o,p);
		q_sort(&a[p+1],&o[p+1],n-p-1);
	}
}

void shade_triangle(triangle *v, color *c){
	
	/*for(int j =0; j<3; j++){
            v->v0[j] -=0.3;
            v->v1[j] +=0.3;
            v->v2[j] -=0.3;
    }*/
	GLfloat a =1;// (GLfloat)rand() / (GLfloat)RAND_MAX;
	glBegin(GL_TRIANGLES);
		glColor4f(c->c0[0],c->c0[1],c->c0[2],a);
        glVertex2i(v->v1[0],v->v1[1]);
		glColor4f(c->c1[0],c->c1[1],c->c1[2],a);
        glVertex2i(v->v1[0],v->v2[1]);
		glColor4f(c->c2[0],c->c2[1],c->c2[2],a);
        glVertex2i(v->v2[0],v->v2[1]);
    glEnd();

}

void render(void)
{
	triangle v;
	color c;
	glEnable( GL_BLEND );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
    glClearColor(1.0,1.0,1.0,0.3);      // Grey ackround color
    glClear(GL_COLOR_BUFFER_BIT);
	int triang_num = sizeof(d.faces) / sizeof(d.faces[0]);

	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0,600.0,600.0,0.0);
	
	// return;
	
	/*GLfloat v0[3] =  {-0.52, 0.34 -0.04};
	GLfloat v1[3] =  {-0.246, -0.65, -0.15};
	GLfloat v2[3] =  {-0.15, -0.50, 0.0};
	GLfloat c0[3] =  {0.742323, 0.229279, 0.356321};
	GLfloat c1[3] =  {0.356321, 0.364300, 0.403574};
	GLfloat c2[3] =  {0.403574, 0.321502, 0.000000};
        
	for(int i = 0; i<3; i++){
		v.v0[i] = v0[i];
		v.v1[i] = v1[i];
		v.v2[i] = v2[i];
		c.c0[i] = c0[i];	
		c.c1[i] = c1[i];
		c.c2[i] = c2[i];
    }*/	
	
	//c.c0[0] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	//c.c0[1] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	//c.c0[2] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	// for(int t=0; t<30; t++){
	// 	for(int i = 0; i<3; i++){
	// 		v.v0[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
	// 		v.v1[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
	// 		v.v2[i] = -1 + (2*((GLfloat)rand()) / (GLfloat)RAND_MAX);
	// 		c.c0[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	// 		c.c1[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
	// 		c.c2[i] = (GLfloat)rand() / (GLfloat)RAND_MAX;
			
	// 	}		
	// 	shade_triangle(&v, &c);
	// }
	int max[2]={0,0};
	int min[2]={600,600};
	for(int t=0; t<triang_num; t++){
		if(t==9999)
			printf("done\n");
		for(int i = 0; i<2; i++){
			v.v0[i] = (GLint)d.verts2d[d.faces[t][0]][i];
			v.v1[i] = (GLint)d.verts2d[d.faces[t][1]][i];
			v.v2[i] = (GLint)d.verts2d[d.faces[t][2]][i];
			
		}
		//printf("t : %d x = %d y = %d\n",t, v.v0[0],v.v0[1]);
		for(int j=0; j<3; j++){	
			c.c0[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			c.c1[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			c.c2[j] = (GLfloat)d.vcolors[d.faces[t][0]][j];
			if(c.c0[j]>0.5)
				printf("t:%d c.c[%d] = %f\n",t,j,c.c0[j]);
			// if(c.c0[j]==0.0 || c.c0[j] == 1.0)
			// 	printf("error");
			// if(c.c1[j]==0.0 || c.c1[j] == 1.0)
			// 	printf("error");
			// if(c.c1[j]==0.0 || c.c2[j] == 1.0)
			// 	printf("error");
			
		}		
		shade_triangle(&v, &c);
		
	}


	glColor4f(0.082353, 0.078431,0.062745,0.1);

	glBegin(GL_TRIANGLES);
        glVertex2i(0,0);
        glVertex2i(500,2);
        glVertex2i(2,500);
    glEnd();

glFlush();
//glutPostRedisplay();
sleep(1);
    
}

int main(int argc, char** argv)
{	
	char const* const fileName = "verts2d.txt"; /* should check that argc > 1 */
	char const* const fileName2 = "vcolors.txt";
	char const* const fileName3 = "faces.txt";
	char const* const fileName4 = "depth.txt";

	int temp_f[10000][3];

    FILE *file , *file2, *file3, *file4;
	// = fopen(fileName, "r"); /* should check the result */
	
    if ((file = fopen(fileName, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file2 = fopen(fileName2, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file3 = fopen(fileName3, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }
	if ((file4 = fopen(fileName4, "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
    }

	read_data_int(4999,2, d.verts2d,file);
	read_data_float(4999, 3, d.vcolors,file2);
	read_data_int(10000, 3, d.faces, file3);
	read_data_float(4999, 1, d.depth, file4);
	//return 0;
	// printf("%d %d\n", d.verts2d[0][0], d.verts2d[0][1]);
	// printf("%d %d\n", d.verts2d[1][0], d.verts2d[1][1]);

	// printf("%d %d %d\n", d.faces[0][0], d.faces[0][1], d.faces[0][2]);
	// printf("%d %d %d\n", d.faces[1][0], d.faces[1][1], d.faces[1][2]);

	int triang_num = sizeof(d.faces) / sizeof(d.faces[0]);
	float *d_Median = (float *) malloc(triang_num * sizeof(float));
	for (int l = 0; l < triang_num; l++) d_Median[l] = 0.0;
	int *d_Order = (int *) malloc(triang_num * sizeof(int)); 
	for (int l = 0; l < triang_num; l++) d_Order[l] = l;

	for(int i=0; i<triang_num;i++){
		for(int j=0; j<3; j++){
			d_Median[i] = d.depth[d.faces[i][j]][0] + d_Median[i]; 
		}
		d_Median[i] = d_Median[i]/3;
	}
	q_sort(d_Median,d_Order,triang_num);

	for(int i=0; i<triang_num;i++)
		for(int j=0; j<3; j++)	
			temp_f[i][j] = d.faces[i][j];

	for(int i=0; i<triang_num;i++)
		for(int j=0; j<3; j++)
			d.faces[i][j] = temp_f[d_Order[i]][j];
			
	// printf("median max = %f \n", d_Median[0]);
	// printf("median min = %f \n", d_Median[triang_num-1]);

    double rgb[3][3];
    float V1[3][2] = {{185.0, 272.0},
                      {355.0, 250.0},
                      {450.0, 300.0}};
    srand(time(0));

    for(int i = 0; i<3; i++){
        for(int j =0; j<3; j++){
            rgb[i][j] = (double)rand() / (double)RAND_MAX;
            //printf(" %f ", rgb[i][j]);
        }
        //printf("\n");
    }
       


    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Fish!");
	
    //init();
    glutDisplayFunc(render);
    glutMainLoop();
	
	fclose(file);
    fclose(file2);
	fclose(file3);
    fclose(file4);
    return 0;
}